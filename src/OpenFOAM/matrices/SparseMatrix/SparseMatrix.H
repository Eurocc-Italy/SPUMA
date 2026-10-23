/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Cineca.
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
    Foam::SparseMatrix

Description
    Abstract interface towards Sparse Matrices.

SourceFiles
    SparseMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SparseMatrix_H
#define Foam_SparseMatrix_H

#include "primitiveFieldsFwd.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

/*---------------------------------------------------------------------------*\
                       Class SparseMatrix Declaration
\*---------------------------------------------------------------------------*/

class SparseMatrix
{
public:

    // Public Types

        //- Enumerated matrix normalisation types
        enum class matrixFormatTypes : char
        {
            DEFAULT_MATRIX_FORMAT,  //!< "default" matrix format (== LDU-COO)
            SELL_MATRIX_FORMAT,     //!< "SELL-C-sigma" matrix format
        };

        //- Names for the matrixFormatTypes
        static const Enum<matrixFormatTypes> matrixFormatTypesNames_;


    // Generated Methods

        //- Default construct
        SparseMatrix() = default;

        //- Copy construct
        SparseMatrix(const SparseMatrix&) = default;

        //- Copy assignment
        SparseMatrix& operator=(const SparseMatrix&) = default;

        virtual ~SparseMatrix() = default;


    // Member Functions

        //- Matrix multiplication with updated interfaces.
        virtual void Amul
        (
            solveScalarField&,
            const tmp<solveScalarField>&,
            const FieldField<Field, scalar>&,
            const lduInterfaceFieldPtrsList&,
            const direction cmpt
        ) const = 0;

        virtual void residual
        (
            solveScalarField& rA,
            const solveScalarField& psi,
            const scalarField& source,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt
        ) const = 0;

        virtual tmp<solveScalarField> residual
        (
            const solveScalarField& psi,
            const scalarField& source,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const direction cmpt
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
