/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Hrvoje Jasak, Wikki Ltd.  All rights reserved
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

    Copyright (C) 2022 Stefano Oliani
    Copyright (C) 2026 Cineca
-------------------------------------------------------------------------------
License
    This file is part of ICSFOAM.

    ICSFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "AMIInterpolation.H"
#include "overlapAMIGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapAMIGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        overlapAMIGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::overlapAMIGAMGInterface::overlapAMIGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface
    (
        index,
        coarseInterfaces
    ),
    fineOverlapAMIInterface_
    (
        refCast<const overlapAMILduInterface>(fineInterface)
    )
{
    // Construct face agglomeration from cell agglomeration
    {
        // From coarse face to cell
        DynamicList<label> dynFaceCells(localRestrictAddressing.size());

        // From face to coarse face
        DynamicList<label> dynFaceRestrictAddressing
        (
            localRestrictAddressing.size()
        );

        Map<label> masterToCoarseFace(localRestrictAddressing.size());

        for (const label curMaster : localRestrictAddressing)
        {
            const auto iter = masterToCoarseFace.cfind(curMaster);

            if (iter.found())
            {
                // Already have coarse face
                dynFaceRestrictAddressing.append(iter.val());
            }
            else
            {
                // New coarse face
                const label coarseI = dynFaceCells.size();
                dynFaceRestrictAddressing.append(coarseI);
                dynFaceCells.append(curMaster);
                masterToCoarseFace.insert(curMaster, coarseI);
            }
        }

        faceCells_.transfer(dynFaceCells);
        faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);
    }


    // On the owner side construct the AMI

    if (fineOverlapAMIInterface_.owner())
    {
        // Construct the neighbour side agglomeration (as the neighbour would
        // do it so it the exact loop above using neighbourRestrictAddressing
        // instead of localRestrictAddressing)

        label nbrCoarseFaceCount = 0;

        labelList nbrFaceRestrictAddressing(poolSwitch(1));
        {
            // From face to coarse face
            DynamicList<label> dynNbrFaceRestrictAddressing
            (
                neighbourRestrictAddressing.size()
            );

            Map<label> masterToCoarseFace(neighbourRestrictAddressing.size());

            for (const label curMaster : neighbourRestrictAddressing)
            {
                const auto iter = masterToCoarseFace.cfind(curMaster);

                if (iter.found())
                {
                    // Already have coarse face
                    dynNbrFaceRestrictAddressing.append(iter.val());
                }
                else
                {
                    // New coarse face
                    const label coarseI = masterToCoarseFace.size();
                    dynNbrFaceRestrictAddressing.append(coarseI);
                    masterToCoarseFace.insert(curMaster, coarseI);
                }
            }

            nbrFaceRestrictAddressing.transfer(dynNbrFaceRestrictAddressing);
            nbrCoarseFaceCount = masterToCoarseFace.size();
        }

        // amiPtr_.reset
        // (
        //     new AMIPatchToPatchInterpolation
        //     (
        //         fineOverlapAMIInterface_.AMI(),
        //         faceRestrictAddressing_,
        //         nbrFaceRestrictAddressing
        //     )
        // );
        const auto expandRestriction =
        [](
            const labelUList& physicalRestriction,
            const label nCopies,
            const label nCoarseFaces
        )
        {
            labelList expanded
            (
                nCopies*physicalRestriction.size(),
                poolSwitch(1)
            );

            for (label copyi = 0; copyi < nCopies; ++copyi)
            {
                forAll(physicalRestriction, facei)
                {
                    expanded[copyi*physicalRestriction.size() + facei] =
                        copyi*nCoarseFaces + physicalRestriction[facei];
                }
            }

            return expanded;
        };

        const auto& fineAMI = fineOverlapAMIInterface_.AMI();

        const labelList expandedSourceRestriction =
            expandRestriction
            (
                faceRestrictAddressing_,
                fineOverlapAMIInterface_.nCopies(),
                faceCells_.size()
            );

        const labelList expandedTargetRestriction =
            expandRestriction
            (
                nbrFaceRestrictAddressing,
                fineOverlapAMIInterface_.neighbPatch().nCopies(),
                nbrCoarseFaceCount
            );

        amiPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                fineAMI,
                expandedSourceRestriction,
                expandedTargetRestriction
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapAMIGAMGInterface::~overlapAMIGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::overlapAMIGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& iF
) const
{
    const overlapAMIGAMGInterface& nbr =
        dynamic_cast<const overlapAMIGAMGInterface&>(neighbPatch());
    const labelUList& nbrFaceCells = nbr.faceCells();

    tmp<labelField> tpnf(new labelField(nbrFaceCells.size()));
    labelField& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = iF[nbrFaceCells[facei]];
    }

    return tpnf;
}


// ************************************************************************* //
