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

Class
    Foam::LduGAMGPreconditioner

Description
    GAMG preconditioner for LduMatrix, based on LduGAMGSolver.

    The preconditioner is a fixed linear operator, as required by Krylov
    solvers such as PBiCGStab: each application performs nVcycles V-cycles
    starting from a zero correction, and the coarsest level is relaxed with
    a fixed number of smoother sweeps (no tolerance). For the same reason
    scaleCorrection and coarsestLevelCorr are not allowed.

    The controls are read from the solver dictionary, e.g.
    \verbatim
    Up
    {
        solver                  PBiCGStab;
        preconditioner          GAMG;
        smoother                twoStageGaussSeidel;
        nVcycles                1;
        nCoarsestSweeps         4;
        nCellsInCoarsestLevel   10;
        ...
    }
    \endverbatim

SourceFiles
    LduGAMGPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LduGAMGPreconditioner_H
#define Foam_LduGAMGPreconditioner_H

#include "LduGAMGSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class LduGAMGPreconditioner Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class LduGAMGPreconditioner
:
    public LduGAMGSolver<Type, DType, LUType>,
    public LduMatrix<Type, DType, LUType>::preconditioner
{
    // Private Data

        //- Number of V-cycles per application
        label nVcycles_;

        //- Number of smoother sweeps on the coarsest level
        label nCoarsestSweeps_;

        //- Coarse-level corrections
        mutable PtrList<Field<Type>> coarseCorrFields_;

        //- Coarse-level sources
        mutable PtrList<Field<Type>> coarseSources_;

        //- Smoothers for all levels
        PtrList<typename LduMatrix<Type, DType, LUType>::smoother> smoothers_;


    // Private Member Functions

        //- No copy construct
        LduGAMGPreconditioner(const LduGAMGPreconditioner&) = delete;

        //- No copy assignment
        void operator=(const LduGAMGPreconditioner&) = delete;


public:

    //- Runtime type information
    TypeName("GAMG");


    // Constructors

        //- Construct from solver and solver controls
        LduGAMGPreconditioner
        (
            const typename LduMatrix<Type, DType, LUType>::solver& sol,
            const dictionary& solverControls
        );


    //- Destructor
    virtual ~LduGAMGPreconditioner() = default;


    // Member Functions

        //- Return wA the preconditioned form of residual rA
        virtual void precondition
        (
            Field<Type>& wA,
            const Field<Type>& rA
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LduGAMGPreconditioner.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
