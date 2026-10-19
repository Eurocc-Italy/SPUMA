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

    foamExecutor exec;

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

        const scalarField sourceCmpt(source.component(cmpt));

        const direction r = row + cmpt;
        const direction c = col + cmpt;
        const scalar* const __restrict__ diagPtr = diag.cbegin();
        const scalar* const __restrict__ sourcePtr = sourceCmpt.cbegin();
        tensor4* __restrict__ blockDiagPtr = blockDiag.begin();
        vector4* __restrict__ blockSourcePtr = blockSource.begin();

        auto LambdaDiag = [=](label celli)
        {
            blockDiagPtr[celli](r, c) = diagPtr[celli];
            blockSourcePtr[celli][r] += sourcePtr[celli];
        };
        exec.parallelFor(LambdaDiag, blockDiag.size());
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        Field<tensor4>& blockUpper = this->upper();

        const scalar* const __restrict__ upperPtr = upper.cbegin();
        tensor4* __restrict__ blockUpperPtr = blockUpper.begin();

        auto LambdaUpper = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockUpperPtr[facei](row + cmpt, col + cmpt) = upperPtr[facei];
            }
        };
        exec.parallelFor(LambdaUpper, upper.size());
    }

    if (matrix.hasLower() || matrix.hasUpper())
    {
        // Symmetric matrices store only upper: const lower() returns it
        const scalarField& lower = std::as_const(matrix).lower();
        Field<tensor4>& blockLower = this->lower();

        const scalar* const __restrict__ lowerPtr = lower.cbegin();
        tensor4* __restrict__ blockLowerPtr = blockLower.begin();

        auto LambdaLower = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockLowerPtr[facei](row + cmpt, col + cmpt) = lowerPtr[facei];
            }
        };
        exec.parallelFor(LambdaLower, lower.size());
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
                new Field<tensor4>(internal.size(), Zero)
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>(internal.size(), Zero)
            );
        }

        const vector* const __restrict__ internalPtr = internal.cbegin();
        const vector* const __restrict__ boundaryPtr = boundary.cbegin();
        tensor4* __restrict__ patchUpperPtr =
            this->interfacesUpper()[patchi].begin();
        tensor4* __restrict__ patchLowerPtr =
            this->interfacesLower()[patchi].begin();

        auto LambdaPatch = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                patchUpperPtr[facei](row + cmpt, col + cmpt) =
                    -boundaryPtr[facei][cmpt];
                patchLowerPtr[facei](row + cmpt, col + cmpt) =
                    internalPtr[facei][cmpt];
            }
        };
        exec.parallelFor(LambdaPatch, internal.size());
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

    foamExecutor exec;

    {
        const scalar* const __restrict__ diagPtr = diag.cbegin();
        const scalar* const __restrict__ sourcePtr = source.cbegin();
        tensor4* __restrict__ blockDiagPtr = blockDiag.begin();
        vector4* __restrict__ blockSourcePtr = blockSource.begin();

        auto LambdaDiag = [=](label celli)
        {
            blockDiagPtr[celli](row, col) = diagPtr[celli];
            blockSourcePtr[celli][row] += sourcePtr[celli];
        };
        exec.parallelFor(LambdaDiag, blockDiag.size());
    }

    if (matrix.hasUpper())
    {
        const scalarField& upper = matrix.upper();
        Field<tensor4>& blockUpper = this->upper();

        const scalar* const __restrict__ upperPtr = upper.cbegin();
        tensor4* __restrict__ blockUpperPtr = blockUpper.begin();

        auto LambdaUpper = [=](label facei)
        {
            blockUpperPtr[facei](row, col) = upperPtr[facei];
        };
        exec.parallelFor(LambdaUpper, upper.size());
    }

    if (matrix.hasLower() || matrix.hasUpper())
    {
        // Symmetric matrices store only upper: const lower() returns it
        const scalarField& lower = std::as_const(matrix).lower();
        Field<tensor4>& blockLower = this->lower();

        const scalar* const __restrict__ lowerPtr = lower.cbegin();
        tensor4* __restrict__ blockLowerPtr = blockLower.begin();

        auto LambdaLower = [=](label facei)
        {
            blockLowerPtr[facei](row, col) = lowerPtr[facei];
        };
        exec.parallelFor(LambdaLower, lower.size());
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
                new Field<tensor4>(internal.size(), Zero)
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>(internal.size(), Zero)
            );
        }

        const scalar* const __restrict__ internalPtr = internal.cbegin();
        const scalar* const __restrict__ boundaryPtr = boundary.cbegin();
        tensor4* __restrict__ patchUpperPtr =
            this->interfacesUpper()[patchi].begin();
        tensor4* __restrict__ patchLowerPtr =
            this->interfacesLower()[patchi].begin();

        auto LambdaPatch = [=](label facei)
        {
            patchUpperPtr[facei](row, col) = -boundaryPtr[facei];
            patchLowerPtr[facei](row, col) = internalPtr[facei];
        };
        exec.parallelFor(LambdaPatch, internal.size());
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

    foamExecutor exec;

    {
        const vector* const __restrict__ diagPtr = matrix.diag().cbegin();
        const vector* const __restrict__ sourcePtr = matrix.source().cbegin();
        tensor4* __restrict__ blockDiagPtr = blockDiag.begin();
        vector4* __restrict__ blockSourcePtr = blockSource.begin();

        auto LambdaDiag = [=](label celli)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockDiagPtr[celli](row + cmpt, col) = diagPtr[celli][cmpt];
                blockSourcePtr[celli][row + cmpt] += sourcePtr[celli][cmpt];
            }
        };
        exec.parallelFor(LambdaDiag, blockDiag.size());
    }

    {
        const vector* const __restrict__ upperPtr = matrix.upper().cbegin();
        const vector* const __restrict__ lowerPtr = matrix.lower().cbegin();
        tensor4* __restrict__ blockUpperPtr = blockUpper.begin();
        tensor4* __restrict__ blockLowerPtr = blockLower.begin();

        auto LambdaOffDiag = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockUpperPtr[facei](row + cmpt, col) = upperPtr[facei][cmpt];
                blockLowerPtr[facei](row + cmpt, col) = lowerPtr[facei][cmpt];
            }
        };
        exec.parallelFor(LambdaOffDiag, blockUpper.size());
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
                    Zero
                )
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>
                (
                    matrix.interfacesLower()[patchi].size(),
                    Zero
                )
            );
        }

        const vector* const __restrict__ srcUpperPtr =
            matrix.interfacesUpper()[patchi].cbegin();
        const vector* const __restrict__ srcLowerPtr =
            matrix.interfacesLower()[patchi].cbegin();
        tensor4* __restrict__ patchUpperPtr =
            this->interfacesUpper()[patchi].begin();
        tensor4* __restrict__ patchLowerPtr =
            this->interfacesLower()[patchi].begin();

        auto LambdaPatch = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                patchUpperPtr[facei](row + cmpt, col) = srcUpperPtr[facei][cmpt];
                patchLowerPtr[facei](row + cmpt, col) = srcLowerPtr[facei][cmpt];
            }
        };
        exec.parallelFor(LambdaPatch, matrix.interfacesUpper()[patchi].size());
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

    foamExecutor exec;

    {
        const vector* const __restrict__ diagPtr = matrix.diag().cbegin();
        const scalar* const __restrict__ sourcePtr = matrix.source().cbegin();
        tensor4* __restrict__ blockDiagPtr = blockDiag.begin();
        vector4* __restrict__ blockSourcePtr = blockSource.begin();

        auto LambdaDiag = [=](label celli)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockDiagPtr[celli](row, col + cmpt) = diagPtr[celli][cmpt];
            }
            blockSourcePtr[celli][row] += sourcePtr[celli];
        };
        exec.parallelFor(LambdaDiag, blockDiag.size());
    }

    {
        const vector* const __restrict__ upperPtr = matrix.upper().cbegin();
        const vector* const __restrict__ lowerPtr = matrix.lower().cbegin();
        tensor4* __restrict__ blockUpperPtr = blockUpper.begin();
        tensor4* __restrict__ blockLowerPtr = blockLower.begin();

        auto LambdaOffDiag = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                blockUpperPtr[facei](row, col + cmpt) = upperPtr[facei][cmpt];
                blockLowerPtr[facei](row, col + cmpt) = lowerPtr[facei][cmpt];
            }
        };
        exec.parallelFor(LambdaOffDiag, blockUpper.size());
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
                    Zero
                )
            );
            this->interfacesLower().set
            (
                patchi,
                new Field<tensor4>
                (
                    matrix.interfacesLower()[patchi].size(),
                    Zero
                )
            );
        }

        const vector* const __restrict__ srcUpperPtr =
            matrix.interfacesUpper()[patchi].cbegin();
        const vector* const __restrict__ srcLowerPtr =
            matrix.interfacesLower()[patchi].cbegin();
        tensor4* __restrict__ patchUpperPtr =
            this->interfacesUpper()[patchi].begin();
        tensor4* __restrict__ patchLowerPtr =
            this->interfacesLower()[patchi].begin();

        auto LambdaPatch = [=](label facei)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                patchUpperPtr[facei](row, col + cmpt) = srcUpperPtr[facei][cmpt];
                patchLowerPtr[facei](row, col + cmpt) = srcLowerPtr[facei][cmpt];
            }
        };
        exec.parallelFor(LambdaPatch, matrix.interfacesUpper()[patchi].size());
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

    Field<vector4> psi(mesh.nCells(), Zero);

    foamExecutor exec;

    vector4* __restrict__ psiPtr = psi.begin();
    {
        const vector* const __restrict__ UPtr = U.primitiveField().cbegin();
        const scalar* const __restrict__ pPtr = p.primitiveField().cbegin();

        auto LambdaPack = [=](label celli)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                psiPtr[celli][UCmpt + cmpt] = UPtr[celli][cmpt];
            }
            psiPtr[celli][pCmpt] = pPtr[celli];
        };
        exec.parallelFor(LambdaPack, psi.size());
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

    {
        vector* __restrict__ UPtr = U.primitiveFieldRef().begin();
        scalar* __restrict__ pPtr = p.primitiveFieldRef().begin();

        auto LambdaUnpack = [=](label celli)
        {
            for (direction cmpt = 0; cmpt < 3; ++cmpt)
            {
                UPtr[celli][cmpt] = psiPtr[celli][UCmpt + cmpt];
            }
            pPtr[celli] = psiPtr[celli][pCmpt];
        };
        exec.parallelFor(LambdaUnpack, psi.size());
    }

    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    return solverPerf;
}

} // End namespace Foam

// ************************************************************************* //