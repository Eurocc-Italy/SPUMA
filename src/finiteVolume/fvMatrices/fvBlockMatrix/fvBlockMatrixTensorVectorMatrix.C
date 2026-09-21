/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2022 Stefano Oliani
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

#include "fvBlockMatrixTensorVectorMatrix.H"
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::fvBlockMatrix<Foam::vector,Foam::tensor>::insertEquation
(
    fvVectorMatrix& matrix
)
{
    // Save a copy for different components
    scalarField& diag = matrix.diag();
    scalarField saveDiag(diag);

    // Add source boundary contribution
    vectorField& source = matrix.source();

    addBoundarySource(matrix, source, false);

    // Get reference to this source field of block system
    vectorField& b = this->source();

    tensorField& blockDiag = this->diag();

    foamExecutor exec;

    for (direction cmptI=0; cmptI<3; cmptI++)
    {
        addBoundaryDiag(matrix, diag, cmptI);

        scalarField sourceCmpt(source.component(cmptI));

        const scalar* const __restrict__ diagPtr = diag.cbegin();
        const scalar* const __restrict__ sourceCmptPtr = sourceCmpt.cbegin();
        tensor* __restrict__ blockDiagPtr = blockDiag.begin();
        vector* __restrict__ bPtr = b.begin();

        auto LambdaDiag = [=](label cellI)
        {
            blockDiagPtr[cellI](cmptI, cmptI) = diagPtr[cellI];
            bPtr[cellI].component(cmptI) += sourceCmptPtr[cellI];
        };
        exec.parallelFor(LambdaDiag, diag.size());

        // Reset diagonal
        diag = saveDiag;
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        tensorField& blockUpper = this->upper();

        const scalar* const __restrict__ upperPtr = upper.cbegin();
        tensor* __restrict__ blockUpperPtr = blockUpper.begin();

        for (direction cmptI=0; cmptI<3; cmptI++)
        {
            auto LambdaUpper = [=](label faceI)
            {
                blockUpperPtr[faceI](cmptI, cmptI) = upperPtr[faceI];
            };
            exec.parallelFor(LambdaUpper, upper.size());
        }
    }

    if (matrix.hasLower())
    {
        const scalarField& lower = matrix.lower();
        tensorField& blockLower = this->lower();

        const scalar* const __restrict__ lowerPtr = lower.cbegin();
        tensor* __restrict__ blockLowerPtr = blockLower.begin();

        for (direction cmptI=0; cmptI<3; cmptI++)
        {
            auto LambdaLower = [=](label faceI)
            {
                blockLowerPtr[faceI](cmptI, cmptI) = lowerPtr[faceI];
            };
            exec.parallelFor(LambdaLower, lower.size());
        }
    }

    const volVectorField& psi = matrix.psi();

    const fvMesh& mesh = psi.mesh();

    this->interfacesUpper().resize(mesh.boundary().size());
    this->interfacesLower().resize(mesh.boundary().size());

    forAll (psi.boundaryField(), patchI)
    {
        const fvPatchField<vector>& pf = psi.boundaryField()[patchI];
        const vectorField& icp = matrix.internalCoeffs()[patchI];
        const vectorField& bcp = matrix.boundaryCoeffs()[patchI];

        this->interfacesUpper().set(patchI, new tensorField(icp.size(), tensor(Zero)));
        this->interfacesLower().set(patchI, new tensorField(icp.size(), tensor(Zero)));

        tensorField& pcoupleUpper = this->interfacesUpper()[patchI];
        tensorField& pcoupleLower = this->interfacesLower()[patchI];

        tensor* __restrict__ pcoupleUpperPtr = pcoupleUpper.begin();
        tensor* __restrict__ pcoupleLowerPtr = pcoupleLower.begin();

        for (direction cmptI=0; cmptI<3; cmptI++)
        {
            scalarField icpCmpt = icp.component(cmptI);
            scalarField bcpCmpt = bcp.component(cmptI);

            const scalar* const __restrict__ icpCmptPtr = icpCmpt.cbegin();
            const scalar* const __restrict__ bcpCmptPtr = bcpCmpt.cbegin();

            auto LambdaPatch = [=](label faceI)
            {
                pcoupleUpperPtr[faceI](cmptI, cmptI) = -bcpCmptPtr[faceI];
                pcoupleLowerPtr[faceI](cmptI, cmptI) = icpCmptPtr[faceI];
            };
            exec.parallelFor(LambdaPatch, pf.size());
        }
    }
}

// ************************************************************************* //
