/*---------------------------------------------------------------------------*\
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

#ifndef Foam_fvBlockMatrixVector4Tensor4Matrix_H
#define Foam_fvBlockMatrixVector4Tensor4Matrix_H

#include "tensor4.H"
#include "vector4.H"

namespace Foam
{

template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    fvVectorMatrix& matrix,
    const label row,
    const label col
);

template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    fvScalarMatrix& matrix,
    const label row,
    const label col
);

template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    const fvBlockMatrix<vector, vector>& matrix,
    const label row,
    const label col
);

template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    const fvBlockMatrix<scalar, vector>& matrix,
    const label row,
    const label col
);

template<>
SolverPerformance<vector4> fvBlockMatrix<vector4, tensor4>::solve
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const label UCmpt,
    GeometricField<scalar, fvPatchField, volMesh>& p,
    const label pCmpt
);

} // End namespace Foam

#endif

// ************************************************************************* //