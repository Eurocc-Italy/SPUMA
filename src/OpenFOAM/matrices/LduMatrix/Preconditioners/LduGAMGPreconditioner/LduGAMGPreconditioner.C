/*---------------------------------------------------------------------------*\
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

#include "LduGAMGPreconditioner.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduGAMGPreconditioner<Type, DType, LUType>::LduGAMGPreconditioner
(
    const typename LduMatrix<Type, DType, LUType>::solver& sol,
    const dictionary& solverControls
)
:
    LduGAMGSolver<Type, DType, LUType>
    (
        sol.fieldName(),
        sol.matrix(),
        solverControls
    ),
    LduMatrix<Type, DType, LUType>::preconditioner(sol),
    nVcycles_(2),
    nCoarsestSweeps_(this->nFinestSweeps_)
{
    this->controlDict_.readIfPresent("nVcycles", nVcycles_);
    this->controlDict_.readIfPresent("nCoarsestSweeps", nCoarsestSweeps_);

    if (nVcycles_ < 1 || nCoarsestSweeps_ < 1)
    {
        FatalIOErrorInFunction(this->controlDict_)
            << "nVcycles and nCoarsestSweeps must be positive"
            << exit(FatalIOError);
    }

    // A fixed linear operator: no residual-dependent scaling and no
    // tolerance-controlled coarsest solver
    if
    (
        this->controlDict_.template getOrDefault<bool>
        (
            "scaleCorrection",
            false
        )
    )
    {
        FatalIOErrorInFunction(this->controlDict_)
            << "scaleCorrection is not allowed in the GAMG preconditioner"
            << exit(FatalIOError);
    }
    this->scaleCorrection_ = false;

    if (this->controlDict_.found("coarsestLevelCorr"))
    {
        FatalIOErrorInFunction(this->controlDict_)
            << "coarsestLevelCorr is not allowed in the GAMG preconditioner:"
            << " use nCoarsestSweeps"
            << exit(FatalIOError);
    }

    const label coarsestLevel = this->matrixLevels_.size() - 1;

    if
    (
        !this->directSolveCoarsest_
     && this->matrixLevels_.set(coarsestLevel)
    )
    {
        dictionary coarsestControls;
        coarsestControls.set("solver", word("SmoothSolver"));
        coarsestControls.set
        (
            "smoother",
            this->controlDict_.template getOrDefault<word>
            (
                "smoother",
                "GaussSeidel"
            )
        );
        // Negative nSweeps: fixed number of sweeps, no residual evaluation
        coarsestControls.set("nSweeps", -nCoarsestSweeps_);

        this->coarsestSolverPtr_ =
            LduMatrix<Type, DType, LUType>::solver::New
            (
                "coarsestLevelCorr",
                this->matrixLevels_[coarsestLevel],
                coarsestControls
            );
    }

    this->initVcycle(coarseCorrFields_, coarseSources_, smoothers_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduGAMGPreconditioner<Type, DType, LUType>::precondition
(
    Field<Type>& wA,
    const Field<Type>& rA
) const
{
    // The finest-level smoother relaxes against the matrix source:
    // replace it with rA while the V-cycles run
    Field<Type>& source =
        const_cast<LduMatrix<Type, DType, LUType>&>(this->matrix_).source();
    const Field<Type> savedSource(source);
    source = rA;

    wA = Zero;

    Field<Type> AwA(wA.size());
    Field<Type> finestCorrection(wA.size());
    Field<Type> finestResidual(rA);

    for (label cycle = 0; cycle < nVcycles_; ++cycle)
    {
        this->Vcycle
        (
            smoothers_,
            wA,
            AwA,
            finestCorrection,
            finestResidual,
            coarseCorrFields_,
            coarseSources_
        );

        if (cycle < nVcycles_ - 1)
        {
            this->matrix_.Amul(AwA, wA);
            finestResidual = rA - AwA;
        }
    }

    source = savedSource;
}


// ************************************************************************* //
