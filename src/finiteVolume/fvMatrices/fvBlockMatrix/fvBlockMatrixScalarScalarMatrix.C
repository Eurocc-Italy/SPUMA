/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2022 Stefano Oliani
    Copyright (C) 2025 Cineca
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


#include "fvBlockMatrixScalarScalarMatrix.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::fvBlockMatrix<Foam::scalar,Foam::scalar>::insertEquation
(
    fvScalarMatrix& matrix
)
{
    scalarField& diag = matrix.diag();
    scalarField& blockDiag = this->diag();

    // Add boundary conditions to diag
    blockDiag = std::move(diag);
    addBoundaryDiag(matrix, blockDiag, 0);

    // Add boundary conditions to source
    scalarField& source = matrix.source();
    addBoundarySource(matrix, source, false);

    foamExecutor exec;

    const scalar* const __restrict__ sourcePtr = source.cbegin();
    scalarField& blockSource = this->source();
    scalar* __restrict__ blockSourcePtr = blockSource.begin();

    auto LambdaDiag = [=](label cellI)
    {
        blockSourcePtr[cellI] += sourcePtr[cellI];
    };
    exec.parallelFor(LambdaDiag, blockDiag.size());

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        scalarField& blockUpper = this->upper();

        blockUpper = std::move(upper);
    }

    if (matrix.hasLower())
    {
        const scalarField& lower = matrix.lower();
        scalarField& blockLower = this->lower();

        blockLower = std::move(lower);
    }

    const volScalarField& psi = matrix.psi();
    const fvMesh& mesh = psi.mesh();

    this->interfacesUpper().resize(mesh.boundary().size());
    this->interfacesLower().resize(mesh.boundary().size());

    forAll (psi.boundaryField(), patchI)
    {
        const scalarField& icp = matrix.internalCoeffs()[patchI];
        const scalarField& bcp = matrix.boundaryCoeffs()[patchI];

        this->interfacesUpper().set(patchI, new scalarField(-bcp));
        this->interfacesLower().set(patchI, new scalarField(icp));
    }
}

// ************************************************************************* //
