/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "SELLMatrix.H"
#include "IOstreams.H"
#include "Switch.H"
#include "objectRegistry.H"
#include "scalarIOField.H"
#include "Time.H"
#include "meshState.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SELLMatrix, 1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::SELLMatrix::SELLMatrix
(
    const lduMatrix& A
)
:
    lduMesh_(A.mesh()),
    lduMatrix_(A)
{
    const SELLAddressing& addr = SELLAddr();

    valuesPtr_ = std::make_unique<scalarField>(addr.numEntries(), Foam::zero{});
    scalar* __restrict__ valuesPtr = (*valuesPtr_).begin();

    const label* const __restrict__ lowerToSELLPtr = addr.lowerToSELL().cbegin();
    const label* const __restrict__ upperToSELLPtr = addr.upperToSELL().cbegin();
    
    const label nInternalFaces = A.lduAddr().upperAddr().size();

    const scalar* const __restrict__ upperPtr = A.upper().cbegin();
    const scalar* const __restrict__ lowerPtr = A.lower().cbegin();

    foamExecutor exec;

    auto LambdaOffDiag = [=](label facei)
    {
        valuesPtr[lowerToSELLPtr[facei]] = lowerPtr[facei];
        valuesPtr[upperToSELLPtr[facei]] = upperPtr[facei];
    };
    exec.parallelFor(LambdaOffDiag, nInternalFaces);
}


Foam::SELLMatrix::SELLMatrix(const SELLMatrix& A)
:
    lduMesh_(A.lduMesh_),
    lduMatrix_(A.lduMatrix_)
{
    if (A.valuesPtr_)
    {
        valuesPtr_ = std::make_unique<scalarField>(*(A.valuesPtr_));
    }
}


Foam::SELLMatrix::SELLMatrix(SELLMatrix&& A)
:
    lduMesh_(A.lduMesh_),
    lduMatrix_(A.lduMatrix_),
    valuesPtr_(std::move(A.valuesPtr_))
{}


Foam::SELLMatrix::SELLMatrix(SELLMatrix& A, bool reuse)
:
    lduMesh_(A.lduMesh_),
    lduMatrix_(A.lduMatrix_)
{
    if (reuse)
    {
        valuesPtr_ = std::move(A.valuesPtr_);
    }
    else
    {
        // Copy assignment
        if (A.valuesPtr_)
        {
            valuesPtr_ = std::make_unique<scalarField>(*(A.valuesPtr_));
        }
    }
}


Foam::SELLMatrix::SELLMatrix
(
    const lduMatrix& A,
    Istream& is
)
:
    lduMesh_(A.mesh()),
    lduMatrix_(A)
{
    bool withValues;

    is >> withValues;

    if (withValues)
    {
        valuesPtr_ = std::make_unique<scalarField>(is);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField& Foam::SELLMatrix::values() const
{
    if (!valuesPtr_)
    {
        FatalErrorInFunction
            << "valuesPtr_ unallocated"
            << abort(FatalError);
    }

    return *valuesPtr_;
}


Foam::scalarField& Foam::SELLMatrix::values(bool init)
{
    if (!valuesPtr_)
    {
        if (init)
        {
            valuesPtr_ =
                std::make_unique<scalarField>(SELLAddr().numEntries(), Foam::zero{});
        }
        else
        {
            valuesPtr_ =
                std::make_unique<scalarField>(SELLAddr().numEntries());
        }
    }

    return *valuesPtr_;
}


Foam::scalarField& Foam::SELLMatrix::values(label size)
{
    if (!valuesPtr_)
    {
        valuesPtr_ = std::make_unique<scalarField>(size, Foam::zero{});
    }

    return *valuesPtr_;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::SELLMatrix::clearOut()
{
    valuesPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const SELLMatrix& mat)
{
    os  << mat.hasValues() << token::SPACE;

    if (mat.hasValues())
    {
        os  << mat.values();
    }

    os.check(FUNCTION_NAME);

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<SELLMatrix>& iproxy
)
{
    const auto& mat = *iproxy;

    os  << "Values:" << Switch::name(mat.hasValues())
        << endl;

    if (mat.hasValues())
    {
        os  << "values:" << mat.values().size() << endl;
    }

    return os;
}


// ************************************************************************* //
