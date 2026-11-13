/*---------------------------------------------------------------------------*\

License
    This file is part of the TurboWG project:
    https://sourceforge.net/projects/turbowg

    The files from the TurboWG project are free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    This source code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this source code.  If not, see <http://www.gnu.org/licenses/>.

    Special notice:
    This file is part of a port of the mixingPlane interface source code from
    the foam-extend project (https://sourceforge.net/projects/foam-extend/) to
    OpenFOAM (https://www.openfoam.com/). Most of the original source code has
    been preserved, including the comments and authorship/IP notices.

Contributors/Copyright:
    Martin Beaudoin, (2023)

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlanePolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::toProfile
(
    const Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    // with communication

    // Algorithm:
    // 1) Master processor holds maps of all zone addressing (data provided)
    // and all remote zone addressing (data required)
    // 2) Each processor will send the locally active data to the master
    // 3) Master assembles all the data
    // 4) Master sends to all processors the data they need to receive
    //    For mixing plane, this is complete profile data.  HJ, 14/Feb/2013
    //
    // Notes:
    // A) If the size of zone addressing is zero, data is not sent
    // B) Communicated data on each processor has the size of live faces
    // C) On processor zero, complete data is available
    // HJ, 4/Jun/2011
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > mixingPlanePolyPatch::toProfile\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    // Prepare profile field to return
    tmp<Field<Type> > tprofileField
    (
        new Field<Type>(nProfileBands(), pTraits<Type>::zero)
    );
    Field<Type>& profileField = tprofileField.ref();

    if (localParallel())
    {
        if (empty())
        {
            // Patch empty, no interpolation
            return tprofileField;
        }

        // Interpolate data to profile.  All data is local
        if (master())
        {
            profileField = patchToPatch().masterToProfile(ff);
        }
        else
        {
            profileField = patchToPatch().slaveToProfile(ff);
        }

        return tprofileField;
    }

#ifdef FOAM_EXTEND_PARALLEL
    if (Pstream::master())
#endif
    {
        // Expand the field to zone size
        Field<Type> expandField(zone().size(), pTraits<Type>::zero);

        // Insert master processor
        const labelList& za = zoneAddressing();

        foamExecutor exec;
        auto expandFieldPtr = expandField.begin();
        const auto zaPtr = za.cbegin();
        const auto ffPtr = ff.cbegin();
        auto Lambda = [=](label i)
        {
            expandFieldPtr[zaPtr[i]] = ffPtr[i];
        };
        exec.parallelFor(Lambda, za.size());

#ifdef FOAM_EXTEND_PARALLEL
        // Master receives and inserts data from all processors for which
        // receiveAddr contains entries
        for (label procI = 1; procI < Pstream::nProcs(); procI++)
        {
            const labelList& curRAddr = receiveAddr()[procI];

            if (!curRAddr.empty())
            {
                Field<Type> receiveBuf(curRAddr.size());

                // Opt: reconsider mode of communication
                UIPstream::read
                (
                    Pstream::commsTypes::blocking,
                    procI,
                    reinterpret_cast<char*>(receiveBuf.begin()),
                    receiveBuf.byteSize()
                );

                // Insert received information
                forAll (curRAddr, i)
                {
                    expandField[curRAddr[i]] = receiveBuf[i];
                }
            }
        }
#endif

        // Interpolate data to profile after parallel expansion
        // Executed on master only to minimise communication
        // HJ, 14/Feb/2013
        if (master())
        {
            profileField = patchToPatch().masterToProfile(expandField);
        }
        else
        {
            profileField = patchToPatch().slaveToProfile(expandField);
        }

#ifdef FOAM_EXTEND_PARALLEL
        // Interpolation to profile complete, send profile to other processors
        for (label procI = 1; procI < Pstream::nProcs(); procI++)
        {
            // Note: profile is sent to all processors, even if patch size
            // is zero.  HJ, 15/Feb/2013

            // Opt: reconsider mode of communication
            OPstream::write
            (
                Pstream::commsTypes::blocking,
                procI,
                reinterpret_cast<const char*>(profileField.begin()),
                profileField.byteSize()
            );
        }
#endif
    }
#ifdef FOAM_EXTEND_PARALLEL
    else
#endif
    {
#ifdef FOAM_EXTEND_PARALLEL
        // Send local data to master and receive remote data
        // If patch is empty, communication is avoided
        // HJ, 4/Jun/2011
        if (!zoneAddressing().empty())
        {
            // Opt: reconsider mode of communication
            OPstream::write
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                reinterpret_cast<const char*>(ff.begin()),
                ff.byteSize()
            );
        }

        // Note: profile is sent to all processors, even if patch size
        // is zero.  HJ, 15/Feb/2013

        // Opt: reconsider mode of communication
        IPstream::read
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo(),
            reinterpret_cast<char*>(profileField.begin()),
            profileField.byteSize()
        );
#endif
    }

    return tprofileField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::mixingPlanePolyPatch::toProfile
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = toProfile(tpf());
    tpf.clear();
    return tint;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::mixingPlanePolyPatch::fromProfile
(
    const Field<Type>& pf
) const
{
    // Create result
    tmp<Field<Type> > tresult
    (
        new Field<Type>(size())
    );

    if (!empty())
    {
        if (master())
        {
            patchToPatch().maskedProfileToMaster
            (
                pf,
                tresult.ref(),
                zoneAddressing()
            );
        }
        else
        {
             patchToPatch().maskedProfileToSlave
             (
                 pf,
                 tresult.ref(),
                 zoneAddressing()
             );
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::mixingPlanePolyPatch::fromProfile
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = fromProfile(tpf());
    tpf.clear();
    return tint;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::interpolate
(
    const Field<Type>& pf
) const
{
    if (localParallel())
    {
        if (empty())
        {
            // Patch empty, no interpolation
            return tmp<Field<Type> >(new Field<Type>());
        }

        if (master())
        {
            return patchToPatch().slaveToMaster(pf);
        }
        else
        {
            return patchToPatch().masterToSlave(pf);
        }
    }
    else
    {
        // Full zone interpolation with parallel communication and masking
        // Note: going from shadow-to-profile-to-local
        const Field<Type> profileField(shadow().toProfile(pf));

        return fromProfile(profileField);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::mixingPlanePolyPatch::interpolate
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = interpolate(tpf());
    tpf.clear();
    return tint;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::mixingPlanePolyPatch::circumferentialAverage
(
    const Field<Type>& pf
) const
{
    if (localParallel())
    {
        if (empty())
        {
            // Patch empty, no interpolation
            return tmp<Field<Type> >(new Field<Type>());
        }

        if (master())
        {
            return patchToPatch().masterToMaster(pf);
        }
        else
        {
            return patchToPatch().slaveToSlave(pf);
        }
    }
    else
    {
        // Full zone interpolation with parallel communication and masking
        // Note: going from local-to-profile-to-local
        return fromProfile(toProfile(pf)());
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::mixingPlanePolyPatch::circumferentialAverage
(
    const tmp<Field<Type> >& tpf
) const
{
    tmp<Field<Type> > tint = circumferentialAverage(tpf());
    tpf.clear();
    return tint;
}


// ************************************************************************* //
