/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
    Copyright (C) 2025 Cineca
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

#include "TtwoStageGaussSeidelSmoother.H"
#include "PrecisionAdaptor.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::TtwoStageGaussSeidelSmoother<Type, DType, LUType>::TtwoStageGaussSeidelSmoother
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix
)
:
    LduMatrix<Type, DType, LUType>::smoother
    (
        fieldName,
        matrix
    ),
    rD_(matrix.diag().size())
{
    const label nCells = matrix.diag().size();
    const DType* const __restrict__ diagPtr = matrix.diag().begin();
    DType* __restrict__ rDPtr = rD_.begin();
    Type localOne(pTraits<Type>::one);

    foamExecutor exec;

    auto LambdaInvDiag = [=](label celli)
    {
        if constexpr(std::is_same<DType, Type>::value)
        {
            rDPtr[celli] = cmptDivide(localOne, diagPtr[celli]);
        }
        else
        {
            rDPtr[celli] = inv(diagPtr[celli]);
        }
    };
    exec.parallelFor(LambdaInvDiag, nCells);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::TtwoStageGaussSeidelSmoother<Type, DType, LUType>::smooth
(
    const word& fieldName_,
    Field<Type>& psi,
    const LduMatrix<Type, DType, LUType>& matrix_,
    const Field<DType>& rD_,
    const label nSweeps
)
{
    Type* __restrict__ psiPtr = psi.begin();
    const Type* const __restrict__ bPtr = matrix_.source().cbegin();

    const label nCells = psi.size();
    const label nInternalFaces = matrix_.upper().size();

    const DType* const __restrict__ diagPtr = matrix_.diag().cbegin();

    const LUType* const __restrict__ upperPtr =
        matrix_.upper().begin();
    const LUType* const __restrict__ lowerPtr =
        matrix_.lower().begin();

    const auto& addr =  matrix_.lduAddr();

    const label* const __restrict__ uPtr =
        addr.upperAddr().begin();
    const label* const __restrict__ lPtr =
        addr.lowerAddr().begin();

    Field<Type> rDr(nCells);
    Type* __restrict__ rDrPtr = rDr.begin();

    Field<Type> gOld(nCells);
    Type* __restrict__ gOldPtr = gOld.begin();

    Field<Type> g(nCells);
    Type* __restrict__ gPtr = g.begin();

    const DType* const __restrict__ rDPtr = rD_.cbegin();

    const scalar omega = 0.9;
    const label nInnerIter = 1;

    foamExecutor exec;

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        // -- Compute new residual vector (scaled by D^-1)

        // --- Calculate A.psi (we use rDr as auxiliary field)
        matrix_.Amul(rDr, psi);

        // --- Calculate rDr = D^-1 * rA
        // --- Initialize g with rDr
        auto LambdarDr = [=](label celli)
        {
            if constexpr(std::is_same<DType, Type>::value)
            {
                const Type rDr = cmptMultiply(rDPtr[celli], (bPtr[celli] - rDrPtr[celli]));
                rDrPtr[celli] = rDr;
                gPtr[celli] = rDr;
                gOldPtr[celli] = rDr;
            }
            else
            {
                const Type rDr = rDPtr[celli] * (bPtr[celli] - rDrPtr[celli]);
                rDrPtr[celli] = rDr;
                gPtr[celli] = rDr;
                gOldPtr[celli] = rDr;
            }
        };
        exec.parallelFor(LambdarDr, nCells);

        // --- Perform local inner (nInnerIter) Jacobi iterations
        for (label j=0; j<nInnerIter; ++j)
        {
            if (j != 0)
            {
                gOld = g;
                g = rDr;
            }

            // --- Multiply g by the matrix omega * D^-1 * L
            if (matrix_.hasLowerCSR())
            {
                const label* const __restrict__ loStartPtr =
                    addr.losortStartAddr().begin();
                const label* const __restrict__ lcsrPtr =
                    addr.lowerCSRAddr().begin();
                const scalar* const __restrict__ lowercsrPtr =
                    matrix_.lowerCSR().begin();

                auto Lambdag = [=](label celli)
                {
                    if constexpr(std::is_same<DType, Type>::value)
                    {
                        const label start = loStartPtr[celli];
                        const label end = loStartPtr[celli + 1];
                        Type val = Zero;

                        for (label i = start; i < end; i++)
                        {
                            const label nbrCell = lcsrPtr[i];
                            val += lowercsrPtr[i] * gOldPtr[nbrCell];
                        }

                        val = cmptMultiply(val, -omega * rDPtr[celli]);
                        gPtr[celli] += val;
                    }
                    else
                    {
                        const label start = loStartPtr[celli];
                        const label end = loStartPtr[celli + 1];
                        Type val = Zero;

                        for (label i = start; i < end; i++)
                        {
                            const label nbrCell = lcsrPtr[i];
                            val += lowercsrPtr[i] * gOldPtr[nbrCell];
                        }

                        val *= -omega * rDPtr[celli];
                        gPtr[celli] += val;
                   }
                };
                exec.parallelFor(Lambdag, nCells);
            }
            else
            {
                auto Lambdag = [=](label facei)
                {
                    if constexpr(std::is_same<DType, Type>::value)
                    {
                        foamAtomic::AtomicAdd
                        (
                            gPtr[uPtr[facei]],
                            cmptMultiply(-omega * rDPtr[uPtr[facei]], lowerPtr[facei] * gOldPtr[lPtr[facei]])
                        );
                    }
                    else
                    {
                        foamAtomic::AtomicAdd
                        (
                            gPtr[uPtr[facei]],
                           -omega * rDPtr[uPtr[facei]] * lowerPtr[facei] * gOldPtr[lPtr[facei]]
                        );
                    }
                };
                exec.parallelFor(Lambdag, nInternalFaces);
            }
        }

        // --- Update solution vector
        psi += g;
    }
}


template<class Type, class DType, class LUType>
void Foam::TtwoStageGaussSeidelSmoother<Type, DType, LUType>::smooth
(
    Field<Type>& psi,
    const label nSweeps
) const
{
    smooth(this->fieldName_, psi, this->matrix_, rD_, nSweeps);
}

// ************************************************************************* //
