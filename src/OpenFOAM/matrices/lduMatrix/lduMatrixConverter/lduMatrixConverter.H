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

Class
    Foam::lduMatrixConverter

Description
    The class is responsible to convert an LDU matrix 
    into a different user-defined one (e.g SELL).

SourceFiles
    lduMatrixConverter.C

\*---------------------------------------------------------------------------*/

#ifndef lduMatrixConverter_H
#define lduMatrixConverter_H

#include "lduMatrix.H"
#include "SparseMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lduMatrixConverter Declaration
\*---------------------------------------------------------------------------*/

class lduMatrixConverter
{

public:

    // Generated Methods

        //- No copy construct
        lduMatrixConverter(const lduMatrixConverter&) = delete;

        //- No copy assignment
        void operator=(const lduMatrixConverter&) = delete;


    //- Constructor
        lduMatrixConverter() = default;


    //- Destructor
        ~lduMatrixConverter() = default;


    // Static Member Functions

        //- Convert an LDU matrix into a different, user-defined one
        static std::shared_ptr<SparseMatrix> convert
        (
            const lduMatrix& mat,
            const SparseMatrix::matrixFormatTypes matrixFormatType
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
