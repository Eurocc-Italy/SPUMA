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
void fvBlockMatrix<vector4, tensor4>::insertMomentumEquation
(
    fvVectorMatrix& matrix
);

template<>
void fvBlockMatrix<vector4, tensor4>::insertPressureEquation
(
    fvScalarMatrix& matrix
);

template<>
void fvBlockMatrix<vector4, tensor4>::insertVelocityPressureCoupling
(
    const fvBlockMatrix<vector, vector>& matrix
);

template<>
void fvBlockMatrix<vector4, tensor4>::insertPressureVelocityCoupling
(
    const fvBlockMatrix<scalar, vector>& matrix
);

} // End namespace Foam

#endif

// ************************************************************************* //