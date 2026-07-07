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

#include "algebraicPairLduGAMGAgglomeration.H"
#include "lduAddressing.H"

const Foam::word& Foam::algebraicPairLduGAMGAgglomeration::typeName()
{
    static const word typeName("algebraicPair");
    return typeName;
}


Foam::algebraicPairLduGAMGAgglomeration::algebraicPairLduGAMGAgglomeration
(
    const lduMesh& mesh,
    const scalarField& faceWeights,
    const dictionary& controlDict
)
:
    LduGAMGAgglomeration(mesh, controlDict),
    mergeLevels_(max(label(1), controlDict.getOrDefault<label>("mergeLevels", 1))),
    forward_(true)
{
    agglomerate(faceWeights);
}


const Foam::word& Foam::algebraicPairLduGAMGAgglomeration::type() const
{
    return typeName();
}


Foam::tmp<Foam::labelField>
Foam::algebraicPairLduGAMGAgglomeration::pairAgglomerate
(
    label& nCoarseCells,
    const lduAddressing& fineMatrixAddressing,
    const scalarField& faceWeights
)
{
    const label nFineCells = fineMatrixAddressing.size();

    const labelUList& upperAddr = fineMatrixAddressing.upperAddr();
    const labelUList& lowerAddr = fineMatrixAddressing.lowerAddr();

    labelList cellFaces(upperAddr.size() + lowerAddr.size());
    labelList cellFaceOffsets(nFineCells + 1);

    {
        labelList nNbrs(nFineCells, Zero);

        forAll(upperAddr, facei)
        {
            nNbrs[upperAddr[facei]]++;
        }

        forAll(lowerAddr, facei)
        {
            nNbrs[lowerAddr[facei]]++;
        }

        cellFaceOffsets[0] = 0;
        forAll(nNbrs, celli)
        {
            cellFaceOffsets[celli + 1] = cellFaceOffsets[celli] + nNbrs[celli];
        }

        nNbrs = 0;

        forAll(upperAddr, facei)
        {
            cellFaces
            [
                cellFaceOffsets[upperAddr[facei]] + nNbrs[upperAddr[facei]]
            ] = facei;

            nNbrs[upperAddr[facei]]++;
        }

        forAll(lowerAddr, facei)
        {
            cellFaces
            [
                cellFaceOffsets[lowerAddr[facei]] + nNbrs[lowerAddr[facei]]
            ] = facei;

            nNbrs[lowerAddr[facei]]++;
        }
    }

    auto tcoarseCellMap = tmp<labelField>::New(nFineCells, -1);
    auto& coarseCellMap = tcoarseCellMap.ref();

    nCoarseCells = 0;

    for (label cellfi = 0; cellfi < nFineCells; ++cellfi)
    {
        const label celli = forward_ ? cellfi : nFineCells - cellfi - 1;

        if (coarseCellMap[celli] < 0)
        {
            label matchFaceNo = -1;
            scalar maxFaceWeight = -GREAT;

            for
            (
                label faceOs = cellFaceOffsets[celli];
                faceOs < cellFaceOffsets[celli + 1];
                ++faceOs
            )
            {
                const label facei = cellFaces[faceOs];

                if
                (
                    coarseCellMap[upperAddr[facei]] < 0
                 && coarseCellMap[lowerAddr[facei]] < 0
                 && faceWeights[facei] > maxFaceWeight
                )
                {
                    matchFaceNo = facei;
                    maxFaceWeight = faceWeights[facei];
                }
            }

            if (matchFaceNo >= 0)
            {
                coarseCellMap[upperAddr[matchFaceNo]] = nCoarseCells;
                coarseCellMap[lowerAddr[matchFaceNo]] = nCoarseCells;
                nCoarseCells++;
            }
            else
            {
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = -GREAT;

                for
                (
                    label faceOs = cellFaceOffsets[celli];
                    faceOs < cellFaceOffsets[celli + 1];
                    ++faceOs
                )
                {
                    const label facei = cellFaces[faceOs];

                    if (faceWeights[facei] > clusterMaxFaceCoeff)
                    {
                        clusterMatchFaceNo = facei;
                        clusterMaxFaceCoeff = faceWeights[facei];
                    }
                }

                if (clusterMatchFaceNo >= 0)
                {
                    coarseCellMap[celli] = max
                    (
                        coarseCellMap[upperAddr[clusterMatchFaceNo]],
                        coarseCellMap[lowerAddr[clusterMatchFaceNo]]
                    );
                }
            }
        }
    }

    for (label cellfi = 0; cellfi < nFineCells; ++cellfi)
    {
        const label celli = forward_ ? cellfi : nFineCells - cellfi - 1;

        if (coarseCellMap[celli] < 0)
        {
            coarseCellMap[celli] = nCoarseCells;
            nCoarseCells++;
        }
    }

    if (!forward_)
    {
        nCoarseCells--;

        forAll(coarseCellMap, celli)
        {
            coarseCellMap[celli] = nCoarseCells - coarseCellMap[celli];
        }

        nCoarseCells++;
    }

    forward_ = !forward_;

    return tcoarseCellMap;
}


void Foam::algebraicPairLduGAMGAgglomeration::agglomerate
(
    const scalarField& startFaceWeights
)
{
    scalarField faceWeights = startFaceWeights;

    label nPairLevels = 0;
    label nCreatedLevels = 0;

    while (nCreatedLevels < maxLevels_ - 1)
    {
        if (!hasMeshLevel(nCreatedLevels))
        {
            FatalErrorInFunction
                << "No mesh at nCreatedLevels:" << nCreatedLevels
                << exit(FatalError);
        }

        const auto& fineMesh = meshLevel(nCreatedLevels);

        label nCoarseCells = -1;

        tmp<labelField> finalAgglomPtr = pairAgglomerate
        (
            nCoarseCells,
            fineMesh.lduAddr(),
            faceWeights
        );

        if
        (
            continueAgglomerating
            (
                nCellsInCoarsestLevel(),
                finalAgglomPtr().size(),
                nCoarseCells,
                fineMesh.comm()
            )
        )
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.set(nCreatedLevels, finalAgglomPtr);
        }
        else
        {
            break;
        }

        agglomerateLduAddressing(nCreatedLevels);

        {
            scalarField aggFaceWeights
            (
                meshLevels_[nCreatedLevels].upperAddr().size(),
                0.0
            );

            restrictFaceField(aggFaceWeights, faceWeights, nCreatedLevels);
            faceWeights = std::move(aggFaceWeights);
        }

        if (nPairLevels % mergeLevels_)
        {
            combineLevels(nCreatedLevels);
        }
        else
        {
            nCreatedLevels++;
        }

        nPairLevels++;
    }

    compactLevels(nCreatedLevels, true);
}


// ************************************************************************* //