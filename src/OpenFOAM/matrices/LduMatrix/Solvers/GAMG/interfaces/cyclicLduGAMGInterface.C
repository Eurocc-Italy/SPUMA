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

#include "cyclicLduGAMGInterface.H"
#include "DynamicList.H"
#include "labelPairHashes.H"
#include "Istream.H"

namespace Foam
{
    defineTypeNameAndDebug(cyclicLduGAMGInterface, 0);
}


Foam::cyclicLduGAMGInterface::cyclicLduGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label,
    const label
)
:
    LduGAMGInterface(index, coarseInterfaces),
    neighbPatchID_
    (
        refCast<const cyclicLduInterface>(fineInterface).neighbPatchID()
    ),
    owner_(refCast<const cyclicLduInterface>(fineInterface).owner()),
    forwardT_(refCast<const cyclicLduInterface>(fineInterface).forwardT()),
    reverseT_(refCast<const cyclicLduInterface>(fineInterface).reverseT())
{
    DynamicList<label> dynamicFaceCells(localRestrictAddressing.size());
    DynamicList<label> dynamicFaceRestrictAddressing
    (
        localRestrictAddressing.size()
    );

    labelPairLookup cellsToCoarseFace(2*localRestrictAddressing.size());

    forAll(localRestrictAddressing, fineFacei)
    {
        labelPair cellPair;

        if (owner())
        {
            cellPair = labelPair
            (
                localRestrictAddressing[fineFacei],
                neighbourRestrictAddressing[fineFacei]
            );
        }
        else
        {
            cellPair = labelPair
            (
                neighbourRestrictAddressing[fineFacei],
                localRestrictAddressing[fineFacei]
            );
        }

        const auto coarseFaceIter = cellsToCoarseFace.cfind(cellPair);

        if (coarseFaceIter.good())
        {
            dynamicFaceRestrictAddressing.append(coarseFaceIter.val());
        }
        else
        {
            const label coarseFacei = dynamicFaceCells.size();
            dynamicFaceRestrictAddressing.append(coarseFacei);
            dynamicFaceCells.append(localRestrictAddressing[fineFacei]);
            cellsToCoarseFace.insert(cellPair, coarseFacei);
        }
    }

    faceCells_.transfer(dynamicFaceCells);
    faceRestrictAddressing_.transfer(dynamicFaceRestrictAddressing);
}


Foam::cyclicLduGAMGInterface::cyclicLduGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const labelUList& faceCells,
    const labelUList& faceRestrictAddressing,
    const label neighbPatchID,
    const bool owner,
    const tensorField& forwardT,
    const tensorField& reverseT
)
:
    LduGAMGInterface(index, coarseInterfaces, faceCells, faceRestrictAddressing),
    neighbPatchID_(neighbPatchID),
    owner_(owner),
    forwardT_(forwardT),
    reverseT_(reverseT)
{}


Foam::cyclicLduGAMGInterface::cyclicLduGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    LduGAMGInterface(index, coarseInterfaces, is),
    neighbPatchID_(readLabel(is)),
    owner_(readBool(is)),
    forwardT_(is),
    reverseT_(is)
{}


Foam::tmp<Foam::labelField>
Foam::cyclicLduGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& internalData
) const
{
    const auto& neighbourPatch =
        refCast<const cyclicLduGAMGInterface>(neighbPatch());
    const labelUList& neighbourFaceCells = neighbourPatch.faceCells();

    auto tpnf = tmp<labelField>::New(size());
    auto& pnf = tpnf.ref();

    forAll(pnf, facei)
    {
        pnf[facei] = internalData[neighbourFaceCells[facei]];
    }

    return tpnf;
}


void Foam::cyclicLduGAMGInterface::write(Ostream& os) const
{
    os  << faceCells_ << token::SPACE << faceRestrictAddressing_
        << token::SPACE << neighbPatchID_
        << token::SPACE << owner_
        << token::SPACE << forwardT_
        << token::SPACE << reverseT_;
}


// ************************************************************************* //