/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
    Copyright (C) 2025-2026 Cineca
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
    Foam::SELLMatrix

Description
    SELLMatrix is a general matrix class designed to optimize performance 
    on modern hardware architectures like GPUs and vector processors (SIMD/SIMT).

    Main features:
    Slices (Chunking): The rows of the sparse matrix are divided into consecutive 
    groups called slices, each containing a fixed number of rows, 
    denoted by C (the chunk size).

    Local Padding: Within each individual slice, the rows are padded with dummy/zero 
    elements so that every row in that specific slice matches the length of the row 
    with the maximum number of non-zeros. This uniformity allows hardware threads 
    to execute instructions synchronously without divergence.

    Data Layout: The non-zero values and their corresponding column indices are stored 
    sequentially in column-major order per slice, paired with an offset array to track 
    where each slice begins.

    Off-diagonals values are stored in the valuesPtr_ array. Diagonal values
    are keep separated and not copied. A reference to lduMatrix is stored.

SourceFiles
    SELLMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SELLMatrix_H
#define Foam_SELLMatrix_H

#include "lduMesh.H"
#include "lduMatrix.H"
#include <functional>  // For reference_wrapper
#include "SELLAddressing.H"
#include "SparseMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class SELLMatrix;


/*---------------------------------------------------------------------------*\
                           Class SELLMatrix Declaration
\*---------------------------------------------------------------------------*/

class SELLMatrix
:
    public SparseMatrix
{
    // Private Data

        //- LDU mesh reference
        std::reference_wrapper<const lduMesh> lduMesh_;

        //- LDU matrix reference
        std::reference_wrapper<const lduMatrix> lduMatrix_;

        //- Off-diagonal coefficients
        std::unique_ptr<scalarField> valuesPtr_;


public:

    // Public Types

    // -----------------------------------------------------------------------

    // Static Data

        // Declare name of the class and its debug switch
        ClassName("SELLMatrix");


    // Constructors

        //- Construct (without coefficients) for an lduMatrix.
        SELLMatrix(const lduMatrix& matrix);

        //- Copy construct
        SELLMatrix(const SELLMatrix&);

        //- Move construct
        SELLMatrix(SELLMatrix&&);

        //- Construct as copy or re-use as specified.
        SELLMatrix(SELLMatrix&, bool reuse);

        //- Construct given an LDU addressed mesh and an Istream
        //- from which the coefficients are read
        SELLMatrix(const lduMatrix& matrix, Istream& is);


    //- Destructor
    virtual ~SELLMatrix() 
    {
        clearOut(); 
    };


    // Private Member Functions

        //- Clear memory
        void clearOut();    


    // Member Functions

    // Addressing

        //- Return the LDU mesh from which the addressing is obtained
        const lduMesh& mesh() const noexcept
        {
            return lduMesh_;
        }

        //- Set the LDU mesh containing the addressing
        void setLduMesh(const lduMesh& m)
        {
            lduMesh_ = m;
        }

        //- Return the SELL addressing
        const SELLAddressing& SELLAddr() const
        {
            return mesh().SELLAddr();
        }

        //- Return the LDU matrix
        const lduMatrix& matrix() const noexcept
        {
            return lduMatrix_;
        }


    // Coefficients

        const scalarField& values() const;

        scalarField& values(bool init = true);

        // Size with externally provided sizes
        // (for constructing with 'fake' mesh in GAMG)

        scalarField& values(label size);


    // Characteristics

        bool hasValues() const noexcept { return bool(valuesPtr_); }

        
        // Operations

            //- Matrix multiplication with updated interfaces.
            void Amul
            (
                solveScalarField&,
                const tmp<solveScalarField>&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;

            //- Matrix transpose multiplication with updated interfaces.
            void Tmul
            (
                solveScalarField&,
                const tmp<solveScalarField>&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&,
                const direction cmpt
            ) const;

            //- Sum the coefficients on each row of the matrix
            void sumA
            (
                solveScalarField&,
                const FieldField<Field, scalar>&,
                const lduInterfaceFieldPtrsList&
            ) const;

            void residual
            (
                solveScalarField& rA,
                const solveScalarField& psi,
                const scalarField& source,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;

            virtual tmp<solveScalarField> residual
            (
                const solveScalarField& psi,
                const scalarField& source,
                const FieldField<Field, scalar>& interfaceBouCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const direction cmpt
            ) const;

            //- Initialise the update of interfaced interfaces
            //- for matrix operations
            void initMatrixInterfaces
            (
                const bool add,
                const FieldField<Field, scalar>& interfaceCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const solveScalarField& psiif,
                solveScalarField& result,
                const direction cmpt
            ) const;

            //- Update interfaced interfaces for matrix operations
            void updateMatrixInterfaces
            (
                const bool add,
                const FieldField<Field, scalar>& interfaceCoeffs,
                const lduInterfaceFieldPtrsList& interfaces,
                const solveScalarField& psiif,
                solveScalarField& result,
                const direction cmpt,
                const label startRequest // starting request (for non-blocking)
            ) const;


    // Ostream Operators

        friend Ostream& operator<<(Ostream&, const SELLMatrix&);
        friend Ostream& operator<<(Ostream&, const InfoProxy<SELLMatrix>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
