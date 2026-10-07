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
#include "vector4ProcessorInterfaceField.H"
#include "vector4CyclicInterfaceField.H"
#include "labelVector.H"

namespace Foam
{

template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    fvVectorMatrix& matrix,
    const label row,
    const label col
)
{
    if (row < 0 || row + 3 > 4 || col < 0 || col + 3 > 4)
    {
        FatalErrorInFunction
            << "Cannot insert a 3x3 equation block at ("
            << row << ',' << col << ") in a tensor4" << abort(FatalError);
    }

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
            blockDiag[celli](row + cmpt, col + cmpt) = diag[celli];
            blockSource[celli][row + cmpt] += source[celli][cmpt];
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
                blockUpper[facei](row + cmpt, col + cmpt) = upper[facei];
            }
        }
    }

    if (matrix.hasLower() || matrix.hasUpper())
    {
        // Symmetric matrices store only upper: const lower() returns it
        const scalarField& lower = std::as_const(matrix).lower();
        Field<tensor4>& blockLower = this->lower();

        forAll(lower, facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockLower[facei](row + cmpt, col + cmpt) = lower[facei];
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
                patchUpper[facei](row + cmpt, col + cmpt) =
                    -boundary[facei][cmpt];
                patchLower[facei](row + cmpt, col + cmpt) =
                    internal[facei][cmpt];
            }
        }
    }
}


template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    fvScalarMatrix& matrix,
    const label row,
    const label col
)
{
    if (row < 0 || row >= 4 || col < 0 || col >= 4)
    {
        FatalErrorInFunction
            << "Cannot insert a 1x1 equation block at ("
            << row << ',' << col << ") in a tensor4" << abort(FatalError);
    }

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
        blockDiag[celli](row, col) = diag[celli];
        blockSource[celli][row] += source[celli];
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        Field<tensor4>& blockUpper = this->upper();

        forAll(upper, facei)
        {
            blockUpper[facei](row, col) = upper[facei];
        }
    }

    if (matrix.hasLower() || matrix.hasUpper())
    {
        // Symmetric matrices store only upper: const lower() returns it
        const scalarField& lower = std::as_const(matrix).lower();
        Field<tensor4>& blockLower = this->lower();

        forAll(lower, facei)
        {
            blockLower[facei](row, col) = lower[facei];
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
            patchUpper[facei](row, col) = -boundary[facei];
            patchLower[facei](row, col) = internal[facei];
        }
    }
}


template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    const fvBlockMatrix<vector, vector>& matrix,
    const label row,
    const label col
)
{
    if (row < 0 || row + 3 > 4 || col < 0 || col >= 4)
    {
        FatalErrorInFunction
            << "Cannot insert a 3x1 equation block at ("
            << row << ',' << col << ") in a tensor4" << abort(FatalError);
    }

    Field<tensor4>& blockDiag = this->diag();
    Field<tensor4>& blockUpper = this->upper();
    Field<tensor4>& blockLower = this->lower();
    Field<vector4>& blockSource = this->source();

    forAll(blockDiag, celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockDiag[celli](row + cmpt, col) = matrix.diag()[celli][cmpt];
            blockSource[celli][row + cmpt] += matrix.source()[celli][cmpt];
        }
    }

    forAll(blockUpper, facei)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockUpper[facei](row + cmpt, col) = matrix.upper()[facei][cmpt];
            blockLower[facei](row + cmpt, col) = matrix.lower()[facei][cmpt];
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
                this->interfacesUpper()[patchi][facei](row + cmpt, col) =
                    matrix.interfacesUpper()[patchi][facei][cmpt];
                this->interfacesLower()[patchi][facei](row + cmpt, col) =
                    matrix.interfacesLower()[patchi][facei][cmpt];
            }
        }
    }
}


template<>
void fvBlockMatrix<vector4, tensor4>::insertEquation
(
    const fvBlockMatrix<scalar, vector>& matrix,
    const label row,
    const label col
)
{
    if (row < 0 || row >= 4 || col < 0 || col + 3 > 4)
    {
        FatalErrorInFunction
            << "Cannot insert a 1x3 equation block at ("
            << row << ',' << col << ") in a tensor4" << abort(FatalError);
    }

    Field<tensor4>& blockDiag = this->diag();
    Field<tensor4>& blockUpper = this->upper();
    Field<tensor4>& blockLower = this->lower();
    Field<vector4>& blockSource = this->source();

    forAll(blockDiag, celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockDiag[celli](row, col + cmpt) = matrix.diag()[celli][cmpt];
        }
        blockSource[celli][row] += matrix.source()[celli];
    }

    forAll(blockUpper, facei)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            blockUpper[facei](row, col + cmpt) = matrix.upper()[facei][cmpt];
            blockLower[facei](row, col + cmpt) = matrix.lower()[facei][cmpt];
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
                this->interfacesUpper()[patchi][facei](row, col + cmpt) =
                    matrix.interfacesUpper()[patchi][facei][cmpt];
                this->interfacesLower()[patchi][facei](row, col + cmpt) =
                    matrix.interfacesLower()[patchi][facei][cmpt];
            }
        }
    }
}


template<>
SolverPerformance<vector4> fvBlockMatrix<vector4, tensor4>::solve
(
    GeometricField<vector, fvPatchField, volMesh>& U,
    const label UCmpt,
    GeometricField<scalar, fvPatchField, volMesh>& p,
    const label pCmpt
)
{
    if (UCmpt < 0 || UCmpt + 3 > 4 || pCmpt < 0 || pCmpt >= 4)
    {
        FatalErrorInFunction
            << "Cannot solve for U at " << UCmpt << " and p at " << pCmpt
            << " in a vector4" << abort(FatalError);
    }

    const fvMesh& mesh = U.mesh();

    Field<vector4> psi(mesh.nCells(), vector4(Zero));

    forAll(psi, celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            psi[celli][UCmpt + cmpt] = U[celli][cmpt];
        }
        psi[celli][pCmpt] = p[celli];
    }

    const lduInterfacePtrsList meshInterfaces = mesh.interfaces();
    PtrList<vector4ProcessorInterfaceField> processorInterfaces
    (
        meshInterfaces.size()
    );
    PtrList<vector4CyclicInterfaceField> cyclicInterfaces
    (
        meshInterfaces.size()
    );
    this->interfaces().resize(meshInterfaces.size());

    forAll(meshInterfaces, patchi)
    {
        if (!meshInterfaces.set(patchi))
        {
            continue;
        }

        if (isA<processorLduInterface>(meshInterfaces[patchi]))
        {
            processorInterfaces.set
            (
                patchi,
                new vector4ProcessorInterfaceField
                (
                    meshInterfaces[patchi],
                    refCast<const processorLduInterfaceField>
                    (
                        U.boundaryField()[patchi]
                    )
                )
            );
            this->interfaces().set(patchi, &processorInterfaces[patchi]);
        }
        else if (isA<cyclicLduInterface>(meshInterfaces[patchi]))
        {
            cyclicInterfaces.set
            (
                patchi,
                new vector4CyclicInterfaceField
                (
                    meshInterfaces[patchi],
                    refCast<const cyclicLduInterfaceField>
                    (
                        U.boundaryField()[patchi]
                    )
                )
            );
            this->interfaces().set(patchi, &cyclicInterfaces[patchi]);
        }
    }

    const word solverName(U.name() + p.name());

    autoPtr<LduMatrix<vector4, tensor4, tensor4>::solver> UpSolver =
        LduMatrix<vector4, tensor4, tensor4>::solver::New
        (
            solverName,
            *this,
            mesh.solverDict(solverName)
        );

    const SolverPerformance<vector4> solverPerf = UpSolver->solve(psi);

    Info<< solverPerf << nl;

    // Store per-field performance so residualControl and solverInfo see U, p
    const vector4& iRes = solverPerf.initialResidual();
    const vector4& fRes = solverPerf.finalResidual();
    const Vector4<label>& nIter = solverPerf.nIterations();

    mesh.data().setSolverPerformance
    (
        U.name(),
        SolverPerformance<vector>
        (
            solverPerf.solverName(),
            U.name(),
            vector(iRes[UCmpt], iRes[UCmpt + 1], iRes[UCmpt + 2]),
            vector(fRes[UCmpt], fRes[UCmpt + 1], fRes[UCmpt + 2]),
            labelVector(nIter[UCmpt], nIter[UCmpt + 1], nIter[UCmpt + 2]),
            solverPerf.converged(),
            solverPerf.singular()
        )
    );

    mesh.data().setSolverPerformance
    (
        p.name(),
        SolverPerformance<scalar>
        (
            solverPerf.solverName(),
            p.name(),
            iRes[pCmpt],
            fRes[pCmpt],
            nIter[pCmpt],
            solverPerf.converged(),
            solverPerf.singular()
        )
    );

    forAll(psi, celli)
    {
        for (direction cmpt = 0; cmpt < 3; ++cmpt)
        {
            U[celli][cmpt] = psi[celli][UCmpt + cmpt];
        }
        p[celli] = psi[celli][pCmpt];
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    return solverPerf;
}

} // End namespace Foam

// ************************************************************************* //