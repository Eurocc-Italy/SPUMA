/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "lduMatrixConverter.H"
#include "SELLMatrix.H"

// * * * * * * * * * * * * * * Static Member Fuctions  * * * * * * * * * * * //

std::shared_ptr<Foam::SparseMatrix> Foam::lduMatrixConverter::convert
(
    const lduMatrix& mat,
    SparseMatrix::matrixFormatTypes matrixFormatType 
)
{
    switch (matrixFormatType)
    {
        case SparseMatrix::matrixFormatTypes::DEFAULT_MATRIX_FORMAT:
        {
            return nullptr;
        }

        case SparseMatrix::matrixFormatTypes::SELL_MATRIX_FORMAT:
        {
            return std::make_shared<SELLMatrix>(mat);
        }

        default:
            FatalErrorInFunction
                << "Unknown matrix format conversion. " 
                << exit(FatalError);
    }
}


// ************************************************************************* //
