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

#include "fvBlockMatrix.H"
#include "fvScalarMatrix.H"

namespace Foam
{

template<>
void fvBlockMatrix<vector4, tensor4>::insertMomentumEquation
(
    fvVectorMatrix& matrix
)
{
    vectorField source(matrix.source());
    forAll(matrix.psi().boundaryField(), patchi)
    {
        if (!matrix.psi().boundaryField()[patchi].coupled())
        {
            addToInternalField
            (
                this->lduAddr().patchAddr(patchi),
                matrix.boundaryCoeffs()[patchi],
                source
            );
        }
    }

    const scalarField internalDiag(matrix.diag());
    Field<tensor4>& blockDiag = this->diag();
    Field<vector4>& blockSource = this->source();

    for (direction cmpt = 0; cmpt < 3; ++cmpt)
    {
        scalarField diag(internalDiag);
        forAll(matrix.internalCoeffs(), patchi)
        {
            addToInternalField
            (
                this->lduAddr().patchAddr(patchi),
                matrix.internalCoeffs()[patchi].component(cmpt),
                diag
            );
        }

        forAll(blockDiag, celli)
        {
            blockDiag[celli](cmpt, cmpt) = diag[celli];
            blockSource[celli][cmpt] += source[celli][cmpt];
        }
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        Field<tensor4>& blockUpper = this->upper();

        forAll(upper, facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockUpper[facei](cmpt, cmpt) = upper[facei];
            }
        }
    }

    if (matrix.hasLower())
    {
        const scalarField& lower = matrix.lower();
        Field<tensor4>& blockLower = this->lower();

        forAll(lower, facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockLower[facei](cmpt, cmpt) = lower[facei];
            }
        }
    }

    const fvMesh& mesh = matrix.psi().mesh();
    this->interfacesUpper().resize(mesh.boundary().size());
    this->interfacesLower().resize(mesh.boundary().size());

    forAll(matrix.psi().boundaryField(), patchi)
    {
        const vectorField& internal = matrix.internalCoeffs()[patchi];
        const vectorField& boundary = matrix.boundaryCoeffs()[patchi];

        if (!this->interfacesUpper().set(patchi))
        {
            this->interfacesUpper().set
            (
                patchi,
                new Field<tensor4>(internal.size(), tensor4(Zero))
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>(internal.size(), tensor4(Zero))
            );
        }

        Field<tensor4>& patchUpper = this->interfacesUpper()[patchi];
        Field<tensor4>& patchLower = this->interfacesLower()[patchi];

        forAll(internal, facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                patchUpper[facei](cmpt, cmpt) = -boundary[facei][cmpt];
                patchLower[facei](cmpt, cmpt) = internal[facei][cmpt];
            }
        }
    }
}


template<>
void fvBlockMatrix<vector4, tensor4>::insertPressureEquation
(
    fvScalarMatrix& matrix
)
{
    scalarField diag(matrix.diag());
    scalarField source(matrix.source());

    forAll(matrix.psi().boundaryField(), patchi)
    {
        addToInternalField
        (
            this->lduAddr().patchAddr(patchi),
            matrix.internalCoeffs()[patchi],
            diag
        );

        if (!matrix.psi().boundaryField()[patchi].coupled())
        {
            addToInternalField
            (
                this->lduAddr().patchAddr(patchi),
                matrix.boundaryCoeffs()[patchi],
                source
            );
        }
    }

    Field<tensor4>& blockDiag = this->diag();
    Field<vector4>& blockSource = this->source();

    forAll(blockDiag, celli)
    {
        blockDiag[celli](3, 3) = diag[celli];
        blockSource[celli][3] += source[celli];
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        Field<tensor4>& blockUpper = this->upper();

        forAll(upper, facei)
        {
            blockUpper[facei](3, 3) = upper[facei];
        }
    }

    if (matrix.hasLower())
    {
        const scalarField& lower = matrix.lower();
        Field<tensor4>& blockLower = this->lower();

        forAll(lower, facei)
        {
            blockLower[facei](3, 3) = lower[facei];
        }
    }

    const fvMesh& mesh = matrix.psi().mesh();
    this->interfacesUpper().resize(mesh.boundary().size());
    this->interfacesLower().resize(mesh.boundary().size());

    forAll(matrix.psi().boundaryField(), patchi)
    {
        const scalarField& internal = matrix.internalCoeffs()[patchi];
        const scalarField& boundary = matrix.boundaryCoeffs()[patchi];

        if (!this->interfacesUpper().set(patchi))
        {
            this->interfacesUpper().set
            (
                patchi,
                new Field<tensor4>(internal.size(), tensor4(Zero))
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>(internal.size(), tensor4(Zero))
            );
        }

        Field<tensor4>& patchUpper = this->interfacesUpper()[patchi];
        Field<tensor4>& patchLower = this->interfacesLower()[patchi];

        forAll(internal, facei)
        {
            patchUpper[facei](3, 3) = -boundary[facei];
            patchLower[facei](3, 3) = internal[facei];
        }
    }
}


template<>
void fvBlockMatrix<vector4, tensor4>::insertVelocityPressureCoupling
(
    const fvBlockMatrix<vector, vector>& matrix
)
{
    Field<tensor4>& blockDiag = this->diag();
    Field<tensor4>& blockUpper = this->upper();
    Field<tensor4>& blockLower = this->lower();
    Field<vector4>& blockSource = this->source();

    forAll(blockDiag, celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockDiag[celli](cmpt, 3) = matrix.diag()[celli][cmpt];
            blockSource[celli][cmpt] += matrix.source()[celli][cmpt];
        }
    }

    forAll(blockUpper, facei)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockUpper[facei](cmpt, 3) = matrix.upper()[facei][cmpt];
            blockLower[facei](cmpt, 3) = matrix.lower()[facei][cmpt];
        }
    }

    this->interfacesUpper().resize(matrix.interfacesUpper().size());
    this->interfacesLower().resize(matrix.interfacesLower().size());

    forAll(matrix.interfacesUpper(), patchi)
    {
        if (!this->interfacesUpper().set(patchi))
        {
            this->interfacesUpper().set
            (
                patchi,
                new Field<tensor4>
                (
                    matrix.interfacesUpper()[patchi].size(),
                    tensor4(Zero)
                )
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>
                (
                    matrix.interfacesLower()[patchi].size(),
                    tensor4(Zero)
                )
            );
        }

        forAll(matrix.interfacesUpper()[patchi], facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                this->interfacesUpper()[patchi][facei](cmpt, 3) =
                    matrix.interfacesUpper()[patchi][facei][cmpt];
                this->interfacesLower()[patchi][facei](cmpt, 3) =
                    matrix.interfacesLower()[patchi][facei][cmpt];
            }
        }
    }
}


template<>
void fvBlockMatrix<vector4, tensor4>::insertPressureVelocityCoupling
(
    const fvBlockMatrix<scalar, vector>& matrix
)
{
    Field<tensor4>& blockDiag = this->diag();
    Field<tensor4>& blockUpper = this->upper();
    Field<tensor4>& blockLower = this->lower();
    Field<vector4>& blockSource = this->source();

    forAll(blockDiag, celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockDiag[celli](3, cmpt) = matrix.diag()[celli][cmpt];
        }
        blockSource[celli][3] += matrix.source()[celli];
    }

    forAll(blockUpper, facei)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockUpper[facei](3, cmpt) = matrix.upper()[facei][cmpt];
            blockLower[facei](3, cmpt) = matrix.lower()[facei][cmpt];
        }
    }

    this->interfacesUpper().resize(matrix.interfacesUpper().size());
    this->interfacesLower().resize(matrix.interfacesLower().size());

    forAll(matrix.interfacesUpper(), patchi)
    {
        if (!this->interfacesUpper().set(patchi))
        {
            this->interfacesUpper().set
            (
                patchi,
                new Field<tensor4>
                (
                    matrix.interfacesUpper()[patchi].size(),
                    tensor4(Zero)
                )
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>
                (
                    matrix.interfacesLower()[patchi].size(),
                    tensor4(Zero)
                )
            );
        }

        forAll(matrix.interfacesUpper()[patchi], facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                this->interfacesUpper()[patchi][facei](3, cmpt) =
                    matrix.interfacesUpper()[patchi][facei][cmpt];
                this->interfacesLower()[patchi][facei](3, cmpt) =
                    matrix.interfacesLower()[patchi][facei][cmpt];
            }
        }
    }
}

} // End namespace Foam

// ************************************************************************* //