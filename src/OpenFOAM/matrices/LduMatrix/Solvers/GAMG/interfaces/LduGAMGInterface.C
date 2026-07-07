/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Cineca
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LduGAMGInterface.H"
#include "cyclicLduGAMGInterface.H"
#include "processorCyclicLduGAMGInterface.H"
#include "processorLduGAMGInterface.H"
#include "Istream.H"

namespace Foam
{
    defineTypeNameAndDebug(LduGAMGInterface, 0);
}


Foam::autoPtr<Foam::LduGAMGInterface> Foam::LduGAMGInterface::New
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
{
    if (fineInterface.type() == processorCyclicLduGAMGInterface::typeName)
    {
        return autoPtr<LduGAMGInterface>
        (
            new processorCyclicLduGAMGInterface
            (
                index,
                coarseInterfaces,
                fineInterface,
                localRestrictAddressing,
                neighbourRestrictAddressing,
                fineLevelIndex,
                coarseComm
            )
        );
    }

    if (isA<processorLduInterface>(fineInterface))
    {
        return autoPtr<LduGAMGInterface>
        (
            new processorLduGAMGInterface
            (
                index,
                coarseInterfaces,
                fineInterface,
                localRestrictAddressing,
                neighbourRestrictAddressing,
                fineLevelIndex,
                coarseComm
            )
        );
    }

    if (isA<cyclicLduInterface>(fineInterface))
    {
        return autoPtr<LduGAMGInterface>
        (
            new cyclicLduGAMGInterface
            (
                index,
                coarseInterfaces,
                fineInterface,
                localRestrictAddressing,
                neighbourRestrictAddressing,
                fineLevelIndex,
                coarseComm
            )
        );
    }

    FatalErrorInFunction
        << "Unsupported LduMatrix GAMG interface type "
        << fineInterface.type()
        << exit(FatalError);

    return nullptr;
}


Foam::autoPtr<Foam::LduGAMGInterface> Foam::LduGAMGInterface::New
(
    const word& coupleType,
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
{
    if (coupleType == processorCyclicLduGAMGInterface::typeName)
    {
        return autoPtr<LduGAMGInterface>
        (
            new processorCyclicLduGAMGInterface(index, coarseInterfaces, is)
        );
    }

    if (coupleType == processorLduGAMGInterface::typeName)
    {
        return autoPtr<LduGAMGInterface>
        (
            new processorLduGAMGInterface(index, coarseInterfaces, is)
        );
    }

    if (coupleType == cyclicLduGAMGInterface::typeName)
    {
        return autoPtr<LduGAMGInterface>
        (
            new cyclicLduGAMGInterface(index, coarseInterfaces, is)
        );
    }

    FatalErrorInFunction
        << "Unsupported LduMatrix GAMG interface type " << coupleType
        << exit(FatalError);

    return nullptr;
}


void Foam::LduGAMGInterface::combine(const LduGAMGInterface& coarseInterface)
{
    const labelList& coarseFaceRestrictAddressing =
        coarseInterface.faceRestrictAddressing_;

    forAll(faceRestrictAddressing_, fineFacei)
    {
        faceRestrictAddressing_[fineFacei] =
            coarseFaceRestrictAddressing[faceRestrictAddressing_[fineFacei]];
    }

    faceCells_ = coarseInterface.faceCells_;
}


Foam::tmp<Foam::labelField> Foam::LduGAMGInterface::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return interfaceInternalField<label>(internalData);
}


Foam::tmp<Foam::labelField> Foam::LduGAMGInterface::interfaceInternalField
(
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    return interfaceInternalField<label>(internalData, faceCells);
}


Foam::tmp<Foam::scalarField> Foam::LduGAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    auto tcoarseCoeffs = tmp<scalarField>::New(size(), Zero);
    auto& coarseCoeffs = tcoarseCoeffs.ref();

    if (fineCoeffs.size() != faceRestrictAddressing_.size())
    {
        FatalErrorInFunction
            << "Size of coefficients " << fineCoeffs.size()
            << " does not correspond to the size of the restriction "
            << faceRestrictAddressing_.size()
            << abort(FatalError);
    }

    if (debug && max(faceRestrictAddressing_) > size())
    {
        FatalErrorInFunction
            << "Face restrict addressing addresses outside of coarse interface"
            << " size. Max addressing:" << max(faceRestrictAddressing_)
            << " coarse size:" << size()
            << abort(FatalError);
    }

    forAll(faceRestrictAddressing_, fineFacei)
    {
        coarseCoeffs[faceRestrictAddressing_[fineFacei]] +=
            fineCoeffs[fineFacei];
    }

    return tcoarseCoeffs;
}


// ************************************************************************* //