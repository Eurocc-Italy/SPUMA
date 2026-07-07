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

#include "processorLduGAMGInterface.H"
#include "DynamicList.H"
#include "labelPairHashes.H"
#include "Istream.H"

namespace Foam
{
    defineTypeNameAndDebug(processorLduGAMGInterface, 0);
}


Foam::processorLduGAMGInterface::processorLduGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label,
    const label coarseComm
)
:
    LduGAMGInterface(index, coarseInterfaces),
    comm_(coarseComm),
    myProcNo_(refCast<const processorLduInterface>(fineInterface).myProcNo()),
    neighbProcNo_
    (
        refCast<const processorLduInterface>(fineInterface).neighbProcNo()
    ),
    forwardT_(refCast<const processorLduInterface>(fineInterface).forwardT()),
    tag_(refCast<const processorLduInterface>(fineInterface).tag())
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

        if (myProcNo() < neighbProcNo())
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


Foam::processorLduGAMGInterface::processorLduGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const labelUList& faceCells,
    const labelUList& faceRestrictAddressing,
    const label comm,
    const label myProcNo,
    const label neighbProcNo,
    const tensorField& forwardT,
    const int tag
)
:
    LduGAMGInterface(index, coarseInterfaces, faceCells, faceRestrictAddressing),
    comm_(comm),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    forwardT_(forwardT),
    tag_(tag)
{}


Foam::processorLduGAMGInterface::processorLduGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    Istream& is
)
:
    LduGAMGInterface(index, coarseInterfaces, is),
    comm_(readLabel(is)),
    myProcNo_(readLabel(is)),
    neighbProcNo_(readLabel(is)),
    forwardT_(is),
    tag_(readLabel(is))
{}


void Foam::processorLduGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& internalData
) const
{
    send(commsType, interfaceInternalField(internalData)());
}


void Foam::processorLduGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& internalData,
    const labelUList& faceCells
) const
{
    send(commsType, interfaceInternalField(internalData, faceCells)());
}


Foam::tmp<Foam::labelField>
Foam::processorLduGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}


void Foam::processorLduGAMGInterface::write(Ostream& os) const
{
    os  << faceCells_ << token::SPACE << faceRestrictAddressing_
        << token::SPACE << comm_
        << token::SPACE << myProcNo_
        << token::SPACE << neighbProcNo_
        << token::SPACE << forwardT_
        << token::SPACE << tag_;
}


// ************************************************************************* //