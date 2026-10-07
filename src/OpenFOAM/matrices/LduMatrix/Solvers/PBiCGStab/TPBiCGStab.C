/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "TPBiCGStab.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::TPBiCGStab<Type, DType, LUType>::TPBiCGStab
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverDict
)
:
    LduMatrix<Type, DType, LUType>::solver
    (
        fieldName,
        matrix,
        solverDict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::SolverPerformance<Type>
Foam::TPBiCGStab<Type, DType, LUType>::solve
(
    Field<Type>& psi
) const
{
    const word preconditionerName(this->controlDict_.getWord("preconditioner"));

    // --- Setup class containing solver performance data
    SolverPerformance<Type> solverPerf
    (
        preconditionerName + typeName,
        this->fieldName_
    );

    const label comm = this->matrix_.mesh().comm();

    label nIter = 0;

    const label nCells = psi.size();

    Type* __restrict__ psiPtr = psi.begin();

    Field<Type> pA(nCells);
    Type* __restrict__ pAPtr = pA.begin();

    Field<Type> yA(nCells);
    Type* __restrict__ yAPtr = yA.begin();

    foamExecutor exec;

    // --- Calculate A.psi
    this->matrix_.Amul(yA, psi);

    // --- Calculate initial residual field
    Field<Type> rA(this->matrix_.source() - yA);
    Type* __restrict__ rAPtr = rA.begin();

    // --- Calculate normalisation factor
    const Type normFactor = this->normFactor(psi, yA, pA);

    if ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Calculate normalised residual norm
    solverPerf.initialResidual() =
        Type(cmptDivide(gSumCmptMag(rA, comm), normFactor));
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Check convergence, solve if not converged
    if
    (
        this->minIter_ > 0
     || !solverPerf.checkConvergence
        (
            this->tolerance_,
            this->relTol_,
            this->log_
        )
    )
    {
        Field<Type> AyA(nCells);
        Type* __restrict__ AyAPtr = AyA.begin();

        Field<Type> sA(nCells);
        Type* __restrict__ sAPtr = sA.begin();

        Field<Type> zA(nCells);
        Type* __restrict__ zAPtr = zA.begin();

        Field<Type> tA(nCells);
        Type* __restrict__ tAPtr = tA.begin();

        // --- Store initial residual
        const Field<Type> rA0(rA);

        // --- Initial values not used
        scalar rA0rA = 0;
        scalar alpha = 0;
        scalar omega = 0;

        // --- Select and construct the preconditioner
        autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        preconPtr = LduMatrix<Type, DType, LUType>::preconditioner::New
        (
            *this,
            this->controlDict_
        );

        // --- Solver iteration
        do
        {
            // --- Store previous rA0rA
            const scalar rA0rAold = rA0rA;

            rA0rA = gSumProd(rA0, rA, comm);

            // --- Test for singularity
            if
            (
                solverPerf.checkSingularity
                (
                    pTraits<Type>::one*mag(rA0rA)
                )
            )
            {
                break;
            }

            // --- Update pA
            if (nIter == 0)
            {
                auto Lambda1 = [=](label cell)
                {
                    pAPtr[cell] = rAPtr[cell];
                };
                exec.parallelFor(Lambda1, nCells);
            }
            else
            {
                // --- Test for singularity
                if
                (
                    solverPerf.checkSingularity
                    (
                        pTraits<Type>::one*mag(omega)
                    )
                )
                {
                    break;
                }

                const scalar beta = (rA0rA/rA0rAold)*(alpha/omega);

                auto Lambda2 = [=](label cell)
                {
                    pAPtr[cell] =
                        rAPtr[cell] + beta*(pAPtr[cell] - omega*AyAPtr[cell]);
                };
                exec.parallelFor(Lambda2, nCells);
            }

            // --- Precondition pA
            preconPtr->precondition(yA, pA);

            // --- Calculate AyA
            this->matrix_.Amul(AyA, yA);

            const scalar rA0AyA = gSumProd(rA0, AyA, comm);

            alpha = rA0rA/rA0AyA;

            // --- Calculate sA
            auto Lambda3 = [=](label cell)
            {
                sAPtr[cell] = rAPtr[cell] - alpha*AyAPtr[cell];
            };
            exec.parallelFor(Lambda3, nCells);

            // --- Test sA for convergence
            solverPerf.finalResidual() =
                Type(cmptDivide(gSumCmptMag(sA, comm), normFactor));

            if
            (
                nIter >= this->minIter_
             && solverPerf.checkConvergence
                (
                    this->tolerance_,
                    this->relTol_,
                    this->log_
                )
            )
            {
                auto Lambda4 = [=](label cell)
                {
                    psiPtr[cell] += alpha*yAPtr[cell];
                };
                exec.parallelFor(Lambda4, nCells);

                ++nIter;

                solverPerf.nIterations() =
                    pTraits<typename pTraits<Type>::labelType>::one*nIter;

                return solverPerf;
            }

            // --- Precondition sA
            preconPtr->precondition(zA, sA);

            // --- Calculate tA
            this->matrix_.Amul(tA, zA);

            const scalar tAtA = gSumProd(tA, tA, comm);

            // --- Calculate omega from tA and sA
            //     (cheaper than using zA with preconditioned tA)
            omega = gSumProd(tA, sA, comm)/tAtA;

            // --- Update solution and residual
            auto Lambda5 = [=](label cell)
            {
                psiPtr[cell] += alpha*yAPtr[cell] + omega*zAPtr[cell];
                rAPtr[cell] = sAPtr[cell] - omega*tAPtr[cell];
            };
            exec.parallelFor(Lambda5, nCells);

            solverPerf.finalResidual() =
                Type(cmptDivide(gSumCmptMag(rA, comm), normFactor));
        } while
        (
            (
                ++nIter < this->maxIter_
             && !solverPerf.checkConvergence
                (
                    this->tolerance_,
                    this->relTol_,
                    this->log_
                )
            )
         || nIter < this->minIter_
        );
    }

    solverPerf.nIterations() =
        pTraits<typename pTraits<Type>::labelType>::one*nIter;

    return solverPerf;
}


// ************************************************************************* //
