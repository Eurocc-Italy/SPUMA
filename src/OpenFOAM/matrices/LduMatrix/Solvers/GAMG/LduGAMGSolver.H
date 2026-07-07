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

Class
    Foam::LduGAMGSolver

Description
    Temporary GAMGC solver scaffold for coupled LduMatrix solvers.

SourceFiles
    LduGAMGSolver.C

\\*---------------------------------------------------------------------------*/

#ifndef Foam_LduGAMGSolver_H
#define Foam_LduGAMGSolver_H

#include "LduGAMGAgglomeration.H"
#include "LduMatrix.H"
#include "scalarMatrices.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\\
                         Class LduGAMGSolver Declaration
\\*---------------------------------------------------------------------------*/


template<class Type, class DType, class LUType>
class LduGAMGSolver
:
    public LduMatrix<Type, DType, LUType>::solver
{

    // Private Data

        //- Number of pre-smoothing sweeps
        label nPreSweeps_;

        //- Lever multiplier for the number of pre-smoothing sweeps
        label preSweepsLevelMultiplier_;

        //- Maximum number of pre-smoothing sweeps
        label maxPreSweeps_;

        //- Number of post-smoothing sweeps
        label nPostSweeps_;

        //- Lever multiplier for the number of post-smoothing sweeps
        label postSweepsLevelMultiplier_;

        //- Maximum number of post-smoothing sweeps
        label maxPostSweeps_;

        //- Number of smoothing sweeps on finest mesh
        label nFinestSweeps_;

        //- Choose if the corrections should be interpolated after injection.
        //  By default corrections are not interpolated.
        bool interpolateCorrection_;

        //- Choose if the corrections should be scaled.
        //  By default corrections for symmetric matrices are scaled
        //  but not for asymmetric matrices.
        bool scaleCorrection_;

        //- Direct or iteratively solve the coarsest level
        bool directSolveCoarsest_;

        //- The native agglomeration hierarchy storage
        autoPtr<LduGAMGAgglomeration> agglomerationPtr_;

        //- The agglomeration hierarchy view
        const LduGAMGAgglomeration& agglomeration_;

        //- Hierarchy of matrix levels
        mutable PtrList<LduMatrix<Type, DType, LUType>> matrixLevels_;

        //- Owning storage for coarse-level interface fields
        PtrList<PtrList<LduInterfaceField<Type>>> primitiveInterfaceLevels_;

        //- Non-owning coarse-level interface views
        PtrList<LduInterfaceFieldPtrsList<Type>> interfaceLevels_;

        //- Coarse-level interface upper coefficients
        PtrList<FieldField<Field, LUType>> interfaceLevelsUpper_;

        //- Coarse-level interface lower coefficients
        PtrList<FieldField<Field, LUType>> interfaceLevelsLower_;

        //- Coarsest-level solver
        autoPtr<typename LduMatrix<Type, DType, LUType>::solver>
            coarsestSolverPtr_;

        //- Communicator over which the coarsest matrix is gathered
        //  for the direct solve
        label coarsestComm_ = -1;

        //- LU-factored gathered coarsest matrix (populated on the master
        //  of coarsestComm_ when directSolveCoarsest is enabled)
        autoPtr<scalarSquareMatrix> coarsestLUPtr_;

        //- Pivot indices of the factored coarsest matrix (master only)
        labelList coarsestPivots_;

        //- Per-rank cell offsets into the gathered coarsest matrix
        labelList coarsestProcOffsets_;

    // Private Member Functions

        //- Return the matrix for a given GAMG level
        const LduMatrix<Type, DType, LUType>& matrixLevel(const label) const;

        //- Return the interface fields for a given GAMG level
        const LduInterfaceFieldPtrsList<Type>& interfaceLevel
        (
            const label
        ) const;

        //- Return the upper interface coefficients for a given GAMG level
        const FieldField<Field, LUType>& interfaceUpperCoeffsLevel
        (
            const label
        ) const;

        //- Return the lower interface coefficients for a given GAMG level
        const FieldField<Field, LUType>& interfaceLowerCoeffsLevel
        (
            const label
        ) const;

        //- Read control parameters from the control dictionary
        void readControls();

        //- Build one coarse matrix level from the corresponding fine level
        void agglomerateMatrix
        (
            const label fineLevelIndex,
            const lduMesh& coarseMesh,
            const lduInterfacePtrsList& coarseMeshInterfaces
        );

        //- Build coarse interface fields and coefficients for one level
        void agglomerateInterfaceCoefficients
        (
            const label fineLevelIndex,
            const lduInterfacePtrsList& coarseMeshInterfaces,
            PtrList<LduInterfaceField<Type>>& coarsePrimInterfaces,
            LduInterfaceFieldPtrsList<Type>& coarseInterfaces,
            FieldField<Field, LUType>& coarseInterfaceUpper,
            FieldField<Field, LUType>& coarseInterfaceLower
        ) const;

        //- Gather matrices (and interface coefficients) from a set of
        //  processors onto the master of the supplied communicator
        void gatherMatrices
        (
            const label destLevel,
            const label comm,

            // Local matrix
            const LduMatrix<Type, DType, LUType>& mat,
            const FieldField<Field, LUType>& interfaceBouCoeffs,
            const FieldField<Field, LUType>& interfaceIntCoeffs,
            const LduInterfaceFieldPtrsList<Type>& interfaces,

            // Remote matrices
            PtrList<LduMatrix<Type, DType, LUType>>& otherMats,
            PtrList<FieldField<Field, LUType>>& otherBouCoeffs,
            PtrList<FieldField<Field, LUType>>& otherIntCoeffs
        ) const;

        //- Assemble the processor-agglomerated matrix on the master
        void procAgglomerateMatrix
        (
            const labelList& procAgglomMap,
            const labelList& agglomProcIDs,
            const label levelI,

            // Resulting matrix
            autoPtr<LduMatrix<Type, DType, LUType>>& allMatrixPtr,
            FieldField<Field, LUType>& allInterfaceBouCoeffs,
            FieldField<Field, LUType>& allInterfaceIntCoeffs,
            PtrList<LduInterfaceField<Type>>& allPrimitiveInterfaces,
            LduInterfaceFieldPtrsList<Type>& allInterfaces
        ) const;

        //- Replace a coarse level with its processor-agglomerated equivalent
        void procAgglomerateMatrix
        (
            const labelList& procAgglomMap,
            const labelList& agglomProcIDs,
            const label levelI
        );

        //- Apply the GAMG correction scaling on one level
        void scale
        (
            Field<Type>& field,
            Field<Type>& Acf,
            const LduMatrix<Type, DType, LUType>& A,
            const Field<Type>& source
        ) const;

        //- Interpolate a correction using the local matrix stencil only
        void interpolate
        (
            Field<Type>& psi,
            Field<Type>& Apsi,
            const LduMatrix<Type, DType, LUType>& A
        ) const;

        //- Interpolate a correction and preserve the coarse correction sum
        void interpolate
        (
            Field<Type>& psi,
            Field<Type>& Apsi,
            const LduMatrix<Type, DType, LUType>& A,
            const labelList& restrictAddressing,
            const Field<Type>& psiC
        ) const;

        //- Initialise temporary fields and smoothers for a V-cycle
        void initVcycle
        (
            PtrList<Field<Type>>& coarseCorrFields,
            PtrList<Field<Type>>& coarseSources,
            PtrList<typename LduMatrix<Type, DType, LUType>::smoother>&
                smoothers
        ) const;

        //- Execute one V-cycle
        void Vcycle
        (
            const PtrList<typename LduMatrix<Type, DType, LUType>::smoother>&
                smoothers,
            Field<Type>& psi,
            Field<Type>& Apsi,
            Field<Type>& finestCorrection,
            Field<Type>& finestResidual,
            PtrList<Field<Type>>& coarseCorrFields,
            PtrList<Field<Type>>& coarseSources
        ) const;

        //- Solve the coarsest correction equation
        void solveCoarsestLevel
        (
            Field<Type>& coarsestCorrField,
            const Field<Type>& coarsestSource
        ) const;

        //- Assemble and LU-factor the coarsest matrix for the direct solve.
        //  Gathers the (possibly distributed) coarsest scalar matrix onto the
        //  master of the coarsest mesh communicator, mirroring the legacy
        //  LUscalarMatrix assembly without depending on it.
        void buildCoarsestDirectMatrix();

        //- Direct (LU) solve of the coarsest correction equation
        void directSolveCoarsestLevel
        (
            Field<Type>& coarsestCorrField,
            const Field<Type>& coarsestSource
        ) const;

        //- No copy construct
        LduGAMGSolver(const LduGAMGSolver&) = delete;

        //- No copy assignment
        void operator=(const LduGAMGSolver&) = delete;


public:

    //- Runtime type information
    TypeName("GAMG");


    // Constructors

        //- Construct from matrix components and solver data dictionary
        LduGAMGSolver
        (
            const word& fieldName,
            const LduMatrix<Type, DType, LUType>& matrix,
            const dictionary& solverControls
        );


    // Destructor

        virtual ~LduGAMGSolver() = default;


    // Member Functions

        //- Solve the matrix with this solver
        virtual SolverPerformance<Type> solve(Field<Type>& psi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LduGAMGSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
