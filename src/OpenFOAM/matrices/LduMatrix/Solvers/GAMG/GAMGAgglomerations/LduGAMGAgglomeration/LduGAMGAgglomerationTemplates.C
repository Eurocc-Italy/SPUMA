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
#include "globalIndex.H"

namespace
{

template<class CoeffType>
Foam::scalarField agglomerationScalarField
(
    const Foam::Field<CoeffType>& coeffs
)
{
    Foam::scalarField result(coeffs.size());

    forAll(coeffs, coeffI)
    {
        result[coeffI] = Foam::mag(coeffs[coeffI]);
    }

    return result;
}


template<class Type, class DType, class LUType>
Foam::scalarField algebraicPairFaceWeights
(
    const Foam::LduMatrix<Type, DType, LUType>& matrix
)
{
    Foam::scalarField upperWeights(agglomerationScalarField(matrix.upper()));

    if (matrix.hasLower())
    {
        return Foam::max
        (
            upperWeights,
            agglomerationScalarField(matrix.lower())
        );
    }

    return upperWeights;
}

} // End anonymous namespace


template<class Type, class DType, class LUType>
Foam::autoPtr<Foam::LduGAMGAgglomeration>
Foam::LduGAMGAgglomeration::New
(
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& controlDict
)
{
    if (!matrix.hasUpper())
    {
        FatalErrorInFunction
            << "LduMatrix GAMG requires an upper coefficient field; "
            << exit(FatalError);
    }

    return LduGAMGAgglomeration::New
    (
        matrix.mesh(),
        algebraicPairFaceWeights(matrix),
        controlDict
    );
}


template<class Type>
void Foam::LduGAMGAgglomeration::restrictField
(
    Field<Type>& coarseField,
    const Field<Type>& fineField,
    const labelList& fineToCoarse
) const
{
    coarseField = Zero;

    forAll(fineField, i)
    {
        coarseField[fineToCoarse[i]] += fineField[i];
    }
}


template<class Type>
void Foam::LduGAMGAgglomeration::restrictField
(
    Field<Type>& coarseField,
    const Field<Type>& fineField,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (!procAgglom && fineField.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << fineField.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    restrictField(coarseField, fineField, fineToCoarse);

    const label coarseLevelIndex = fineLevelIndex + 1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        if (UPstream::myProcNo(coarseComm) == (procIDs.empty() ? 0 : procIDs[0]))
        {
            coarseField.resize_nocopy(offsets.back());
        }

        globalIndex::gather
        (
            offsets,
            coarseComm,
            procIDs,
            coarseField,
            coarseField,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );
    }
}


template<class Type>
void Foam::LduGAMGAgglomeration::restrictFaceField
(
    Field<Type>& coarseField,
    const Field<Type>& fineField,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (fineField.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << fineField.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    coarseField = Zero;

    forAll(fineToCoarse, fineFacei)
    {
        const label coarseFacei = fineToCoarse[fineFacei];

        if (coarseFacei >= 0)
        {
            coarseField[coarseFacei] += fineField[fineFacei];
        }
    }
}


template<class Type>
void Foam::LduGAMGAgglomeration::prolongField
(
    Field<Type>& fineField,
    const Field<Type>& coarseField,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    const label coarseLevelIndex = fineLevelIndex + 1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        const label localSize = nCells_[fineLevelIndex];

        Field<Type> allCf(localSize);
        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            coarseField,
            allCf,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );

        forAll(fineToCoarse, i)
        {
            fineField[i] = allCf[fineToCoarse[i]];
        }
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            fineField[i] = coarseField[fineToCoarse[i]];
        }
    }
}


template<class Type>
const Foam::Field<Type>& Foam::LduGAMGAgglomeration::prolongField
(
    Field<Type>& fineField,
    Field<Type>& work,
    const Field<Type>& coarseField,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    const label coarseLevelIndex = fineLevelIndex + 1;

    if (hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        const label localSize = nCells_[fineLevelIndex];
        work.resize_nocopy(localSize);

        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            coarseField,
            work,
            UPstream::msgType(),
            UPstream::commsTypes::nonBlocking
        );

        forAll(fineToCoarse, i)
        {
            fineField[i] = work[fineToCoarse[i]];
        }
        return work;
    }
    else
    {
        forAll(fineToCoarse, i)
        {
            fineField[i] = coarseField[fineToCoarse[i]];
        }
        return coarseField;
    }
}


// ************************************************************************* //