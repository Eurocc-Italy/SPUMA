/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
    Copyright (C) 2025-2026 Cineca
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

/*  Added */
#ifdef SPUMA_SELL_SVE
#include <arm_sve.h>
#include <cstdint>
#include <cstddef>
#include "debug.H"
#include "registerSwitch.H"

// Loop unrolling factor, chosen at compile time:
// only the corresponding kernel is compiled.
// Set by Make/options from the environment variable
// LOOP_UNROLLING_FACTOR_G (1, 2 or 4; default 2).
#ifndef LOOP_UNROLLING_FACTOR_G
#define LOOP_UNROLLING_FACTOR_G 2
#endif

// Only the intrinsics kernels are compiled with the SVE target attribute,
// the rest of the code is compiled for the default target.
#define SELL_SVE_TARGET __attribute__((target("+sve")))

namespace
{
    //- Optimisation switch sellVector (default 1): 0 selects the SPUMA loop in SELLMatrixATmul.C, 1 selects the SVE kernels (C > 1).
    int sellVectorSwitch(Foam::debug::optimisationSwitch("sellVector", 1));

    //- True if the SVE kernels are used (sellVector != 0 and C > 1).
    //  The first call reports C and G.
    bool sellVectorKernels(const Foam::label C)
    {
        static const bool on = (sellVectorSwitch != 0);
        static bool reported = false;
        if (on && C > 1 && !reported)
        {
            Foam::Info<< "[SELLMatrixATmul.C] SVE kernels, C = " << C
                << ", G = " << LOOP_UNROLLING_FACTOR_G << Foam::endl;
            reported = true;
        }
        return on && C > 1;
    }

    // The kernels take 32-bit labels and double precision
    static_assert
    (
        sizeof(Foam::label) == sizeof(int),
        "the SELL SVE kernels take 32-bit indices (WM_LABEL_SIZE=32)"
    );
    static_assert
    (
        sizeof(Foam::scalar) == sizeof(double)
     && sizeof(Foam::solveScalar) == sizeof(double),
        "the SELL SVE kernels are double precision values (WM_PRECISION_OPTION=DP)"
    );
}
registerOptSwitch("sellVector", int, sellVectorSwitch);
#endif
/* /Added */

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
    const label C = SELLAddressing::C;

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
    /*  Added:
        SELL Amul vectorized in SVE with fixed G factor (LOOP_UNROLLING_FACTOR_G) and C > 1, otherwise SPUMA loop.
    */
    #ifdef SPUMA_SELL_SVE
    if (sellVectorKernels(C))
    {
        const uint32_t* const __restrict__ boffPtr = addr.boff().data();
        const label nSlices = addr.numSlices();

        // The kernel of G = LOOP_UNROLLING_FACTOR_G (the only one compiled)
        #if LOOP_UNROLLING_FACTOR_G == 1
        [&]() SELL_SVE_TARGET
        {
            // mv_sellcsigma_u32_fused_unroll_1
            const int n = nCells;
            const int nslices = nSlices;
            const int* const slice_ptr = sliceBeginsPtr;
            const int* const slice_width = sliceWidthsPtr;
            const double* const diag_perm = diagPtr;
            const double* const sell_values = valuesPtr;
            const uint32_t* const sell_boff = boffPtr;
            const double* const x_perm = psiPtr;
            double* const y_perm = ApsiPtr;

            const uint64_t vl = svcntd();   // doubles per register
            for (int s = 0; s < nslices; ++s)
            {
                const int      r0   = s * C;
                const uint64_t rows = (uint64_t)((r0 + C > n) ? (n - r0) : C);
                const int      w    = slice_width[s];
                const size_t   bp   = (size_t)slice_ptr[s];
                double        *ys   = &y_perm[r0];

                for (uint64_t rl = 0; rl < rows; rl += vl)
                {
                    svbool_t pg = svwhilelt_b64_u64(rl, rows);   // tail via predicate
                    svfloat64_t d  = svld1_f64(pg, &diag_perm[r0 + rl]);
                    svfloat64_t xv = svld1_f64(pg, &x_perm[r0 + rl]);
                    svst1_f64(pg, &ys[rl], svmul_f64_x(pg, d, xv));
                }

                for (int slot = 0; slot < w; ++slot)
                {
                    const size_t p = bp + (size_t)slot * (size_t)C;
                    for (uint64_t rl = 0; rl < rows; rl += vl)
                    {
                        svbool_t pg = svwhilelt_b64_u64(rl, rows);   // tail via predicate

                        svfloat64_t vals = svld1_f64(pg, &sell_values[p + rl]);
                        svuint64_t  boff = svld1uw_u64(pg, &sell_boff[p + rl]);   // u32 zero-extended to 64-bit lanes
                        svbool_t    mask = svcmpne_n_u64(pg, boff, UINT32_MAX);   // padding sentinel; pg also folds in the tail
                        svfloat64_t g    = svld1_gather_u64offset_f64(mask, x_perm, boff);   // gather with 64-bit byte offsets

                        svfloat64_t y = svld1_f64(pg, &ys[rl]);
                        y = svmla_f64_m(mask, y, vals, g);   // masked-off lanes keep y
                        svst1_f64(pg, &ys[rl], y);
                    }
                }
            }
        }();
        #elif LOOP_UNROLLING_FACTOR_G == 2
        [&]() SELL_SVE_TARGET
        {
            // mv_sellcsigma_u32_fused_unroll_2
            const int n = nCells;
            const int nslices = nSlices;
            const int* const slice_ptr = sliceBeginsPtr;
            const int* const slice_width = sliceWidthsPtr;
            const double* const diag_perm = diagPtr;
            const double* const sell_values = valuesPtr;
            const uint32_t* const sell_boff = boffPtr;
            const double* const x_perm = psiPtr;
            double* const y_perm = ApsiPtr;

            const uint64_t vl = svcntd();
            for (int s = 0; s < nslices; ++s)
            {
                const int      r0   = s * C;
                const uint64_t rows = (uint64_t)((r0 + C > n) ? (n - r0) : C);
                const int      w    = slice_width[s];
                const size_t   bp   = (size_t)slice_ptr[s];
                double        *ys   = &y_perm[r0];
                const double  *ds   = &diag_perm[r0];
                const double  *xs   = &x_perm[r0];

                for (uint64_t rl = 0; rl < rows; rl += 2 * vl)
                {
                    svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                    svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);

                    svfloat64_t d0 = svld1_f64(p0, &ds[rl]);
                    svfloat64_t d1 = svld1_f64(p1, &ds[rl + vl]);
                    svfloat64_t x0 = svld1_f64(p0, &xs[rl]);
                    svfloat64_t x1 = svld1_f64(p1, &xs[rl + vl]);
                    svst1_f64(p0, &ys[rl],      svmul_f64_x(p0, d0, x0));
                    svst1_f64(p1, &ys[rl + vl], svmul_f64_x(p1, d1, x1));
                }

                for (int slot = 0; slot < w; ++slot)
                {
                    const double   *vs = &sell_values[bp + (size_t)slot * (size_t)C];
                    const uint32_t *bs = &sell_boff[bp + (size_t)slot * (size_t)C];
                    for (uint64_t rl = 0; rl < rows; rl += 2 * vl)
                    {
                        svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                        svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);

                        svfloat64_t v0 = svld1_f64(p0, &vs[rl]);
                        svfloat64_t v1 = svld1_f64(p1, &vs[rl + vl]);
                        svuint64_t  b0 = svld1uw_u64(p0, &bs[rl]);
                        svuint64_t  b1 = svld1uw_u64(p1, &bs[rl + vl]);
                        svbool_t    m0 = svcmpne_n_u64(p0, b0, UINT32_MAX);
                        svbool_t    m1 = svcmpne_n_u64(p1, b1, UINT32_MAX);
                        svfloat64_t g0 = svld1_gather_u64offset_f64(m0, x_perm, b0);
                        svfloat64_t g1 = svld1_gather_u64offset_f64(m1, x_perm, b1);

                        svfloat64_t y0 = svld1_f64(p0, &ys[rl]);
                        svfloat64_t y1 = svld1_f64(p1, &ys[rl + vl]);
                        svst1_f64(p0, &ys[rl],      svmla_f64_m(m0, y0, v0, g0));
                        svst1_f64(p1, &ys[rl + vl], svmla_f64_m(m1, y1, v1, g1));
                    }
                }
            }
        }();
        #elif LOOP_UNROLLING_FACTOR_G == 4
        [&]() SELL_SVE_TARGET
        {
            // mv_sellcsigma_u32_fused_unroll_4
            const int n = nCells;
            const int nslices = nSlices;
            const int* const slice_ptr = sliceBeginsPtr;
            const int* const slice_width = sliceWidthsPtr;
            const double* const diag_perm = diagPtr;
            const double* const sell_values = valuesPtr;
            const uint32_t* const sell_boff = boffPtr;
            const double* const x_perm = psiPtr;
            double* const y_perm = ApsiPtr;

            const uint64_t vl = svcntd();
            for (int s = 0; s < nslices; ++s)
            {
                const int      r0   = s * C;
                const uint64_t rows = (uint64_t)((r0 + C > n) ? (n - r0) : C);
                const int      w    = slice_width[s];
                const size_t   bp   = (size_t)slice_ptr[s];
                double        *ys   = &y_perm[r0];
                const double  *ds   = &diag_perm[r0];
                const double  *xs   = &x_perm[r0];

                for (uint64_t rl = 0; rl < rows; rl += 4 * vl)
                {
                    svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                    svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);
                    svbool_t p2 = svwhilelt_b64_u64(rl + 2 * vl, rows);
                    svbool_t p3 = svwhilelt_b64_u64(rl + 3 * vl, rows);

                    svfloat64_t d0 = svld1_f64(p0, &ds[rl]);
                    svfloat64_t d1 = svld1_f64(p1, &ds[rl + vl]);
                    svfloat64_t d2 = svld1_f64(p2, &ds[rl + 2 * vl]);
                    svfloat64_t d3 = svld1_f64(p3, &ds[rl + 3 * vl]);
                    svfloat64_t x0 = svld1_f64(p0, &xs[rl]);
                    svfloat64_t x1 = svld1_f64(p1, &xs[rl + vl]);
                    svfloat64_t x2 = svld1_f64(p2, &xs[rl + 2 * vl]);
                    svfloat64_t x3 = svld1_f64(p3, &xs[rl + 3 * vl]);
                    svst1_f64(p0, &ys[rl],          svmul_f64_x(p0, d0, x0));
                    svst1_f64(p1, &ys[rl + vl],     svmul_f64_x(p1, d1, x1));
                    svst1_f64(p2, &ys[rl + 2 * vl], svmul_f64_x(p2, d2, x2));
                    svst1_f64(p3, &ys[rl + 3 * vl], svmul_f64_x(p3, d3, x3));
                }

                for (int slot = 0; slot < w; ++slot)
                {
                    const double   *vs = &sell_values[bp + (size_t)slot * (size_t)C];
                    const uint32_t *bs = &sell_boff[bp + (size_t)slot * (size_t)C];
                    for (uint64_t rl = 0; rl < rows; rl += 4 * vl)
                    {
                        svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                        svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);
                        svbool_t p2 = svwhilelt_b64_u64(rl + 2 * vl, rows);
                        svbool_t p3 = svwhilelt_b64_u64(rl + 3 * vl, rows);

                        svfloat64_t v0 = svld1_f64(p0, &vs[rl]);
                        svfloat64_t v1 = svld1_f64(p1, &vs[rl + vl]);
                        svfloat64_t v2 = svld1_f64(p2, &vs[rl + 2 * vl]);
                        svfloat64_t v3 = svld1_f64(p3, &vs[rl + 3 * vl]);
                        svuint64_t  b0 = svld1uw_u64(p0, &bs[rl]);
                        svuint64_t  b1 = svld1uw_u64(p1, &bs[rl + vl]);
                        svuint64_t  b2 = svld1uw_u64(p2, &bs[rl + 2 * vl]);
                        svuint64_t  b3 = svld1uw_u64(p3, &bs[rl + 3 * vl]);
                        svbool_t    m0 = svcmpne_n_u64(p0, b0, UINT32_MAX);
                        svbool_t    m1 = svcmpne_n_u64(p1, b1, UINT32_MAX);
                        svbool_t    m2 = svcmpne_n_u64(p2, b2, UINT32_MAX);
                        svbool_t    m3 = svcmpne_n_u64(p3, b3, UINT32_MAX);
                        svfloat64_t g0 = svld1_gather_u64offset_f64(m0, x_perm, b0);
                        svfloat64_t g1 = svld1_gather_u64offset_f64(m1, x_perm, b1);
                        svfloat64_t g2 = svld1_gather_u64offset_f64(m2, x_perm, b2);
                        svfloat64_t g3 = svld1_gather_u64offset_f64(m3, x_perm, b3);

                        svfloat64_t y0 = svld1_f64(p0, &ys[rl]);
                        svfloat64_t y1 = svld1_f64(p1, &ys[rl + vl]);
                        svfloat64_t y2 = svld1_f64(p2, &ys[rl + 2 * vl]);
                        svfloat64_t y3 = svld1_f64(p3, &ys[rl + 3 * vl]);
                        svst1_f64(p0, &ys[rl],          svmla_f64_m(m0, y0, v0, g0));
                        svst1_f64(p1, &ys[rl + vl],     svmla_f64_m(m1, y1, v1, g1));
                        svst1_f64(p2, &ys[rl + 2 * vl], svmla_f64_m(m2, y2, v2, g2));
                        svst1_f64(p3, &ys[rl + 3 * vl], svmla_f64_m(m3, y3, v3, g3));
                    }
                }
            }
        }();
        #else
        #error "LOOP_UNROLLING_FACTOR_G must be 1, 2 or 4"
        #endif
    }
    else
    #endif
    /* /Added */
    {
        exec.parallelFor(LambdaAmul, nCells);
    }

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
    const label C = SELLAddressing::C;

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
    const label C = SELLAddressing::C;

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
    /*  Added: the same residual with SVE intrinsics. */
    #ifdef SPUMA_SELL_SVE
    if (sellVectorKernels(C))
    {
        const uint32_t* const __restrict__ boffPtr = addr.boff().data();
        const label nSlices = addr.numSlices();

        // The kernel of G = LOOP_UNROLLING_FACTOR_G (the only one compiled)
        #if LOOP_UNROLLING_FACTOR_G == 1
        [&]() SELL_SVE_TARGET
        {
            // mv_sellcsigma_u32_fused_res_unroll_1
            const int n = nCells;
            const int nslices = nSlices;
            const int* const slice_ptr = sliceBeginsPtr;
            const int* const slice_width = sliceWidthsPtr;
            const double* const diag_perm = diagPtr;
            const double* const sell_values = valuesPtr;
            const uint32_t* const sell_boff = boffPtr;
            const double* const x_perm = psiPtr;
            const double* const b_perm = sourcePtr;
            double* const r_perm = rAPtr;

            const uint64_t vl = svcntd();   // doubles per register
            for (int s = 0; s < nslices; ++s)
            {
                const int      r0   = s * C;
                const uint64_t rows = (uint64_t)((r0 + C > n) ? (n - r0) : C);
                const int      w    = slice_width[s];
                const size_t   bp   = (size_t)slice_ptr[s];
                double        *ys   = &r_perm[r0];   // this slice's r

                for (uint64_t rl = 0; rl < rows; rl += vl)
                {
                    svbool_t pg = svwhilelt_b64_u64(rl, rows);   // tail via predicate
                    svfloat64_t d  = svld1_f64(pg, &diag_perm[r0 + rl]);
                    svfloat64_t c  = svld1_f64(pg, &b_perm[r0 + rl]);
                    svfloat64_t xv = svld1_f64(pg, &x_perm[r0 + rl]);
                    svst1_f64(pg, &ys[rl], svmls_f64_x(pg, c, d, xv));
                }

                for (int slot = 0; slot < w; ++slot)
                {
                    const size_t p = bp + (size_t)slot * (size_t)C;
                    for (uint64_t rl = 0; rl < rows; rl += vl)
                    {
                        svbool_t pg = svwhilelt_b64_u64(rl, rows);   // tail via predicate

                        svfloat64_t vals = svld1_f64(pg, &sell_values[p + rl]);
                        svuint64_t  boff = svld1uw_u64(pg, &sell_boff[p + rl]);   // u32 zero-extended to 64-bit lanes
                        svbool_t    mask = svcmpne_n_u64(pg, boff, UINT32_MAX);   // padding sentinel; pg also folds in the tail
                        svfloat64_t g    = svld1_gather_u64offset_f64(mask, x_perm, boff);   // gather with 64-bit byte offsets

                        svfloat64_t y = svld1_f64(pg, &ys[rl]);   // current r (b - d*x minus the slots before this one)
                        y = svmls_f64_m(mask, y, vals, g);        // _m merges: padding lanes (mask off) keep r, so their garbage vals and unloaded g never reach it
                        svst1_f64(pg, &ys[rl], y);                // pg, not mask: padding lanes store the r they loaded, unchanged
                    }
                }
            }
        }();
        #elif LOOP_UNROLLING_FACTOR_G == 2
        [&]() SELL_SVE_TARGET
        {
            // mv_sellcsigma_u32_fused_res_unroll_2
            const int n = nCells;
            const int nslices = nSlices;
            const int* const slice_ptr = sliceBeginsPtr;
            const int* const slice_width = sliceWidthsPtr;
            const double* const diag_perm = diagPtr;
            const double* const sell_values = valuesPtr;
            const uint32_t* const sell_boff = boffPtr;
            const double* const x_perm = psiPtr;
            const double* const b_perm = sourcePtr;
            double* const r_perm = rAPtr;

            const uint64_t vl = svcntd();
            for (int s = 0; s < nslices; ++s)
            {
                const int      r0   = s * C;
                const uint64_t rows = (uint64_t)((r0 + C > n) ? (n - r0) : C);
                const int      w    = slice_width[s];
                const size_t   bp   = (size_t)slice_ptr[s];
                double        *ys   = &r_perm[r0];   // this slice's r
                const double  *ds   = &diag_perm[r0];
                const double  *bsrc = &b_perm[r0];
                const double  *xs   = &x_perm[r0];

                for (uint64_t rl = 0; rl < rows; rl += 2 * vl)
                {
                    svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                    svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);

                    svfloat64_t d0 = svld1_f64(p0, &ds[rl]);
                    svfloat64_t c0 = svld1_f64(p0, &bsrc[rl]);
                    svfloat64_t d1 = svld1_f64(p1, &ds[rl + vl]);
                    svfloat64_t c1 = svld1_f64(p1, &bsrc[rl + vl]);
                    svfloat64_t x0 = svld1_f64(p0, &xs[rl]);
                    svfloat64_t x1 = svld1_f64(p1, &xs[rl + vl]);
                    svst1_f64(p0, &ys[rl],      svmls_f64_x(p0, c0, d0, x0));
                    svst1_f64(p1, &ys[rl + vl], svmls_f64_x(p1, c1, d1, x1));
                }

                for (int slot = 0; slot < w; ++slot)
                {
                    const double   *vs = &sell_values[bp + (size_t)slot * (size_t)C];
                    const uint32_t *bs = &sell_boff[bp + (size_t)slot * (size_t)C];
                    for (uint64_t rl = 0; rl < rows; rl += 2 * vl)
                    {
                        svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                        svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);

                        svfloat64_t v0 = svld1_f64(p0, &vs[rl]);
                        svfloat64_t v1 = svld1_f64(p1, &vs[rl + vl]);
                        svuint64_t  b0 = svld1uw_u64(p0, &bs[rl]);
                        svuint64_t  b1 = svld1uw_u64(p1, &bs[rl + vl]);
                        svbool_t    m0 = svcmpne_n_u64(p0, b0, UINT32_MAX);
                        svbool_t    m1 = svcmpne_n_u64(p1, b1, UINT32_MAX);
                        svfloat64_t g0 = svld1_gather_u64offset_f64(m0, x_perm, b0);
                        svfloat64_t g1 = svld1_gather_u64offset_f64(m1, x_perm, b1);

                        svfloat64_t y0 = svld1_f64(p0, &ys[rl]);
                        svfloat64_t y1 = svld1_f64(p1, &ys[rl + vl]);
                        svst1_f64(p0, &ys[rl],      svmls_f64_m(m0, y0, v0, g0));
                        svst1_f64(p1, &ys[rl + vl], svmls_f64_m(m1, y1, v1, g1));
                    }
                }
            }
        }();
        #elif LOOP_UNROLLING_FACTOR_G == 4
        [&]() SELL_SVE_TARGET
        {
            // mv_sellcsigma_u32_fused_res_unroll_4
            const int n = nCells;
            const int nslices = nSlices;
            const int* const slice_ptr = sliceBeginsPtr;
            const int* const slice_width = sliceWidthsPtr;
            const double* const diag_perm = diagPtr;
            const double* const sell_values = valuesPtr;
            const uint32_t* const sell_boff = boffPtr;
            const double* const x_perm = psiPtr;
            const double* const b_perm = sourcePtr;
            double* const r_perm = rAPtr;

            const uint64_t vl = svcntd();
            for (int s = 0; s < nslices; ++s)
            {
                const int      r0   = s * C;
                const uint64_t rows = (uint64_t)((r0 + C > n) ? (n - r0) : C);
                const int      w    = slice_width[s];
                const size_t   bp   = (size_t)slice_ptr[s];
                double        *ys   = &r_perm[r0];   // this slice's r
                const double  *ds   = &diag_perm[r0];
                const double  *bsrc = &b_perm[r0];
                const double  *xs   = &x_perm[r0];

                for (uint64_t rl = 0; rl < rows; rl += 4 * vl)
                {
                    svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                    svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);
                    svbool_t p2 = svwhilelt_b64_u64(rl + 2 * vl, rows);
                    svbool_t p3 = svwhilelt_b64_u64(rl + 3 * vl, rows);

                    svfloat64_t d0 = svld1_f64(p0, &ds[rl]);
                    svfloat64_t c0 = svld1_f64(p0, &bsrc[rl]);
                    svfloat64_t d1 = svld1_f64(p1, &ds[rl + vl]);
                    svfloat64_t c1 = svld1_f64(p1, &bsrc[rl + vl]);
                    svfloat64_t d2 = svld1_f64(p2, &ds[rl + 2 * vl]);
                    svfloat64_t c2 = svld1_f64(p2, &bsrc[rl + 2 * vl]);
                    svfloat64_t d3 = svld1_f64(p3, &ds[rl + 3 * vl]);
                    svfloat64_t c3 = svld1_f64(p3, &bsrc[rl + 3 * vl]);
                    svfloat64_t x0 = svld1_f64(p0, &xs[rl]);
                    svfloat64_t x1 = svld1_f64(p1, &xs[rl + vl]);
                    svfloat64_t x2 = svld1_f64(p2, &xs[rl + 2 * vl]);
                    svfloat64_t x3 = svld1_f64(p3, &xs[rl + 3 * vl]);
                    svst1_f64(p0, &ys[rl],          svmls_f64_x(p0, c0, d0, x0));
                    svst1_f64(p1, &ys[rl + vl],     svmls_f64_x(p1, c1, d1, x1));
                    svst1_f64(p2, &ys[rl + 2 * vl], svmls_f64_x(p2, c2, d2, x2));
                    svst1_f64(p3, &ys[rl + 3 * vl], svmls_f64_x(p3, c3, d3, x3));
                }

                for (int slot = 0; slot < w; ++slot)
                {
                    const double   *vs = &sell_values[bp + (size_t)slot * (size_t)C];
                    const uint32_t *bs = &sell_boff[bp + (size_t)slot * (size_t)C];
                    for (uint64_t rl = 0; rl < rows; rl += 4 * vl)
                    {
                        svbool_t p0 = svwhilelt_b64_u64(rl, rows);
                        svbool_t p1 = svwhilelt_b64_u64(rl + vl, rows);
                        svbool_t p2 = svwhilelt_b64_u64(rl + 2 * vl, rows);
                        svbool_t p3 = svwhilelt_b64_u64(rl + 3 * vl, rows);

                        svfloat64_t v0 = svld1_f64(p0, &vs[rl]);
                        svfloat64_t v1 = svld1_f64(p1, &vs[rl + vl]);
                        svfloat64_t v2 = svld1_f64(p2, &vs[rl + 2 * vl]);
                        svfloat64_t v3 = svld1_f64(p3, &vs[rl + 3 * vl]);
                        svuint64_t  b0 = svld1uw_u64(p0, &bs[rl]);
                        svuint64_t  b1 = svld1uw_u64(p1, &bs[rl + vl]);
                        svuint64_t  b2 = svld1uw_u64(p2, &bs[rl + 2 * vl]);
                        svuint64_t  b3 = svld1uw_u64(p3, &bs[rl + 3 * vl]);
                        svbool_t    m0 = svcmpne_n_u64(p0, b0, UINT32_MAX);
                        svbool_t    m1 = svcmpne_n_u64(p1, b1, UINT32_MAX);
                        svbool_t    m2 = svcmpne_n_u64(p2, b2, UINT32_MAX);
                        svbool_t    m3 = svcmpne_n_u64(p3, b3, UINT32_MAX);
                        svfloat64_t g0 = svld1_gather_u64offset_f64(m0, x_perm, b0);
                        svfloat64_t g1 = svld1_gather_u64offset_f64(m1, x_perm, b1);
                        svfloat64_t g2 = svld1_gather_u64offset_f64(m2, x_perm, b2);
                        svfloat64_t g3 = svld1_gather_u64offset_f64(m3, x_perm, b3);

                        svfloat64_t y0 = svld1_f64(p0, &ys[rl]);
                        svfloat64_t y1 = svld1_f64(p1, &ys[rl + vl]);
                        svfloat64_t y2 = svld1_f64(p2, &ys[rl + 2 * vl]);
                        svfloat64_t y3 = svld1_f64(p3, &ys[rl + 3 * vl]);
                        svst1_f64(p0, &ys[rl],          svmls_f64_m(m0, y0, v0, g0));
                        svst1_f64(p1, &ys[rl + vl],     svmls_f64_m(m1, y1, v1, g1));
                        svst1_f64(p2, &ys[rl + 2 * vl], svmls_f64_m(m2, y2, v2, g2));
                        svst1_f64(p3, &ys[rl + 3 * vl], svmls_f64_m(m3, y3, v3, g3));
                    }
                }
            }
        }();
        #else
        #error "LOOP_UNROLLING_FACTOR_G must be 1, 2 or 4"
        #endif
    }
    else
    #endif
    /* /Added */
    {
        exec.parallelFor(LambdaResidual, nCells);
    }

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
