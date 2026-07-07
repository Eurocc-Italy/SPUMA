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

#include "LduGAMGAgglomeration.H"
#include "algebraicPairLduGAMGAgglomeration.H"
#include "dictionary.H"

namespace Foam
{

LduGAMGAgglomeration::LduGAMGAgglomeration
(
    const lduMesh& mesh,
    const dictionary& controlDict
)
:
    maxLevels_(50),
    nCellsInCoarsestLevel_
    (
        controlDict.getOrDefault<label>("nCellsInCoarsestLevel", 10)
    ),
    mesh_(mesh),
    meshInterfaces_(mesh.interfaces()),
    nCells_(maxLevels_),
    restrictAddressing_(maxLevels_),
    nFaces_(maxLevels_),
    faceRestrictAddressing_(maxLevels_),
    faceFlipMap_(maxLevels_),
    nPatchFaces_(maxLevels_),
    patchFaceRestrictAddressing_(maxLevels_),
    meshLevels_(maxLevels_)
{
    nCellsInCoarsestLevel_ =
        max(1, min(mesh.lduAddr().size()/2, nCellsInCoarsestLevel_));

    reduce(nCellsInCoarsestLevel_, minOp<label>());

    // Processor agglomeration controls (masterCoarsest strategy only)
    doProcessorAgglomerate_ = false;
    nProcessorsPerMaster_ = 0;

    const word procAgglomeratorType
    (
        controlDict.getOrDefault<word>("processorAgglomerator", "none")
    );

    if (procAgglomeratorType == "masterCoarsest")
    {
        doProcessorAgglomerate_ = true;

        nProcessorsPerMaster_ =
            controlDict.getOrDefault<label>
            (
                "nProcessorsPerMaster",
                0,
                keyType::LITERAL
            );

        const auto* ePtr =
            controlDict.findEntry("nMasters", keyType::LITERAL);

        if (ePtr)
        {
            if (nProcessorsPerMaster_ > 0)
            {
                FatalIOErrorInFunction(controlDict)
                    << "Cannot specify both \"nMasters\" and"
                    << " \"nProcessorsPerMaster\"" << exit(FatalIOError);
            }

            const label nMasters(readLabel(ePtr->stream()));

            if (nMasters <= 0)
            {
                FatalIOErrorInFunction(controlDict)
                    << "Illegal value \"nMasters\" "
                    << nMasters << exit(FatalIOError);
            }

            nProcessorsPerMaster_ =
                (UPstream::nProcs(mesh.comm()) + nMasters - 1) / nMasters;
        }

        if (nProcessorsPerMaster_ < 0)
        {
            FatalIOErrorInFunction(controlDict)
                << "Illegal value \"nProcessorsPerMaster\" "
                << nProcessorsPerMaster_ << exit(FatalIOError);
        }
    }
    else if (procAgglomeratorType != "none")
    {
        FatalErrorInFunction
            << "Unsupported LduMatrix GAMG processor agglomerator '"
            << procAgglomeratorType << "'. Only 'masterCoarsest' is"
            << " implemented in the native path" << exit(FatalError);
    }
}


autoPtr<LduGAMGAgglomeration> LduGAMGAgglomeration::New
(
    const lduMesh& mesh,
    const scalarField& faceWeights,
    const dictionary& controlDict
)
{
    const word agglomeratorType
    (
        controlDict.getOrDefault<word>
        (
            "agglomerator",
            algebraicPairLduGAMGAgglomeration::typeName()
        )
    );

    if (agglomeratorType != algebraicPairLduGAMGAgglomeration::typeName())
    {
        FatalErrorInFunction
            << "Unsupported LduMatrix GAMG agglomerator '" << agglomeratorType
            << "'. Only '" << algebraicPairLduGAMGAgglomeration::typeName()
            << "' is implemented in the native serial path"
            << exit(FatalError);
    }

    return autoPtr<LduGAMGAgglomeration>
    (
        new algebraicPairLduGAMGAgglomeration(mesh, faceWeights, controlDict)
    );
}


void LduGAMGAgglomeration::compactLevels
(
    const label nCreatedLevels,
    const bool doProcessorAgglomerate
)
{
    nCells_.setSize(nCreatedLevels, 0);
    restrictAddressing_.setSize(nCreatedLevels);
    nFaces_.setSize(nCreatedLevels, 0);
    faceRestrictAddressing_.setSize(nCreatedLevels);
    faceFlipMap_.setSize(nCreatedLevels);
    nPatchFaces_.setSize(nCreatedLevels);
    patchFaceRestrictAddressing_.setSize(nCreatedLevels);
    meshLevels_.setSize(nCreatedLevels);

    // Have procCommunicator_ always, even if not procAgglomerating
    procCommunicator_.setSize(nCreatedLevels + 1, -1);

    if (doProcessorAgglomerate && doProcessorAgglomerate_)
    {
        procAgglomMap_.setSize(nCreatedLevels);
        agglomProcIDs_.setSize(nCreatedLevels);
        procAgglomCommunicator_.setSize(nCreatedLevels);
        procCellOffsets_.setSize(nCreatedLevels);
        procFaceMap_.setSize(nCreatedLevels);
        procBoundaryMap_.setSize(nCreatedLevels);
        procBoundaryFaceMap_.setSize(nCreatedLevels);

        procAgglomerate();
    }
}


bool LduGAMGAgglomeration::continueAgglomerating
(
    const label nCellsInCoarsestLevel,
    const label nFineCells,
    const label nCoarseCells,
    const label comm
) const
{
    const label nTotalCoarseCells =
        returnReduce(nCoarseCells, sumOp<label>(), UPstream::msgType(), comm);

    if (nTotalCoarseCells < Pstream::nProcs(comm)*nCellsInCoarsestLevel)
    {
        return false;
    }

    const label nTotalFineCells =
        returnReduce(nFineCells, sumOp<label>(), UPstream::msgType(), comm);

    return nTotalCoarseCells < nTotalFineCells;
}


const lduMesh& LduGAMGAgglomeration::meshLevel(const label leveli) const
{
    return leveli ? meshLevels_[leveli - 1] : mesh_;
}


bool LduGAMGAgglomeration::hasMeshLevel(const label leveli) const
{
    return !leveli || meshLevels_.set(leveli - 1);
}


const lduInterfacePtrsList& LduGAMGAgglomeration::interfaceLevel
(
    const label leveli
) const
{
    return leveli ? meshLevels_[leveli - 1].rawInterfaces() : meshInterfaces_;
}

} // End namespace Foam


// ************************************************************************* //