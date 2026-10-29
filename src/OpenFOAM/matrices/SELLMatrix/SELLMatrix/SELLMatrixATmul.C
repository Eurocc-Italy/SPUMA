/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
    Copyright (C) 2026 Cineca
    Copyright (C) 2026 Gabriele Ceccolini - Unibo
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

Description
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "SELLMatrix.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::SELLMatrix::Amul
(
    solveScalarField& Apsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ ApsiPtr = Apsi.begin();

    const solveScalarField& psi = tpsi();
    const solveScalar* const __restrict__ psiPtr = psi.begin();

    const auto& addr = SELLAddr();
    
    const scalar* const __restrict__ diagPtr = matrix().diag().begin();
    const scalar* const __restrict__ valuesPtr = values().begin();

    const label nCells = psi.size();
    const label C = addr.C();

    const label* const __restrict__ sliceWidthsPtr = addr.sliceWidths().cbegin();
    const label* const __restrict__ sliceBeginsPtr = addr.sliceBegins().cbegin();
    const label* const __restrict__ colsPtr = addr.cols().cbegin();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    matrix().initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt
    );

    foamExecutor exec;

    auto LambdaAmul = [=](label rowi)
    {
        const label slicei = rowi / C;    // blockIdx.x
        const label sliceRowi = rowi % C; // threadIdx.x
        
        const label width = sliceWidthsPtr[slicei];
        const label sliceBegin = sliceBeginsPtr[slicei];

        solveScalar sum = diagPtr[rowi] * psiPtr[rowi];

        for (label slot = 0; slot < width; ++slot)
        {
            const label slotRowi = sliceBegin + slot * C + sliceRowi;
            const label col = colsPtr[slotRowi];

            // Check if it's actual data or zero-padding
            if (col != INT32_MAX)
            {
                sum += valuesPtr[slotRowi] * psiPtr[col];
            }
        }

        ApsiPtr[rowi] = sum; 
    };
    exec.parallelFor(LambdaAmul, nCells);

    //Update interface interfaces
    matrix().updateMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        Apsi,
        cmpt,
        startRequest
    );

    tpsi.clear();
}


void Foam::SELLMatrix::Tmul
(
    solveScalarField& Tpsi,
    const tmp<solveScalarField>& tpsi,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    NotImplemented;
}


void Foam::SELLMatrix::sumA
(
    solveScalarField& sumA,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    solveScalar* __restrict__ sumAPtr = sumA.begin();

    const auto& addr = SELLAddr();

    const scalar* const __restrict__ diagPtr = matrix().diag().begin();
    const scalar* const __restrict__ valuesPtr = values().begin();

    const label nCells = matrix().diag().size();
    const label C = addr.C();

    const label* const __restrict__ sliceWidthsPtr = addr.sliceWidths().cbegin();
    const label* const __restrict__ sliceBeginsPtr = addr.sliceBegins().cbegin();
    const label* const __restrict__ colsPtr = addr.cols().cbegin();

    foamExecutor exec;

    auto LambdaSumA = [=](label rowi)
    {
        const label slicei = rowi / C;    // blockIdx.x
        const label sliceRowi = rowi % C; // threadIdx.x

        const label width = sliceWidthsPtr[slicei];
        const label sliceBegin = sliceBeginsPtr[slicei];

        solveScalar sumA = diagPtr[rowi];

        for (label slot = 0; slot < width; ++slot)
        {
            const label slotRowi = sliceBegin + slot * C + sliceRowi;
            const label col = colsPtr[slotRowi];

            // Check if it's actual data or zero-padding
            if (col != INT32_MAX)
            {
                sumA += valuesPtr[slotRowi];
            }
        }

        sumAPtr[rowi] = sumA;
    };
    exec.parallelFor(LambdaSumA, nCells);

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces, patchi)
    {
        if (interfaces.set(patchi))
        {
            const labelUList& pa = matrix().lduAddr().patchAddr(patchi);
            const scalarField& pCoeffs = interfaceBouCoeffs[patchi];

            const auto pap = pa.cbegin();
            const auto pCoeffsp = pCoeffs.cbegin();
            auto Lambda = [=](label face)
            {
                foamAtomic::AtomicAdd(sumAPtr[pap[face]], -pCoeffsp[face]);
            };
            exec.parallelFor(Lambda, pa.size());
        }
    }
}


void Foam::SELLMatrix::residual
(
    solveScalarField& rA,
    const solveScalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    solveScalar* __restrict__ rAPtr = rA.begin();

    const solveScalar* const __restrict__ psiPtr = psi.begin();
    const solveScalar* const __restrict__ sourcePtr = source.cbegin();

    const auto& addr = SELLAddr();

    const scalar* const __restrict__ diagPtr = matrix().diag().begin();
    const scalar* const __restrict__ valuesPtr = values().begin();

    const label nCells = psi.size();
    const label C = addr.C();

    const label* const __restrict__ sliceWidthsPtr = addr.sliceWidths().cbegin();
    const label* const __restrict__ sliceBeginsPtr = addr.sliceBegins().cbegin();
    const label* const __restrict__ colsPtr = addr.cols().cbegin();

    const label startRequest = UPstream::nRequests();

    // Initialise the update of interfaced interfaces
    matrix().initMatrixInterfaces
    (
        true,
        interfaceBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt
    );

    foamExecutor exec;

    auto LambdaResidual = [=](label rowi)
    {
        const label slicei = rowi / C;    // blockIdx.x
        const label sliceRowi = rowi % C; // threadIdx.x

        const label width = sliceWidthsPtr[slicei];
        const label sliceBegin = sliceBeginsPtr[slicei];

        solveScalar res = sourcePtr[rowi] - diagPtr[rowi] * psiPtr[rowi];

        for (label slot = 0; slot < width; ++slot)
        {
            const label slotRowi = sliceBegin + slot * C + sliceRowi;
            const label col = colsPtr[slotRowi];

            // Check if it's actual data or zero-padding
            if (col != INT32_MAX)
            {
                res -= valuesPtr[slotRowi] * psiPtr[col];
            }
        }

        rAPtr[rowi] = res;
    };
    exec.parallelFor(LambdaResidual, nCells);

    //Update interface interfaces
    matrix().updateMatrixInterfaces
    (
        false,
        interfaceBouCoeffs,
        interfaces,
        psi,
        rA,
        cmpt,
        startRequest
    );
}


Foam::tmp<Foam::Field<Foam::solveScalar>> Foam::SELLMatrix::residual
(
    const solveScalarField& psi,
    const scalarField& source,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    auto trA = tmp<solveScalarField>::New(psi.size());
    residual(trA.ref(), psi, source, interfaceBouCoeffs, interfaces, cmpt);
    return trA;
}


// ************************************************************************* //
