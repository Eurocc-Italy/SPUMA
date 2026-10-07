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

#include "vector4CyclicInterfaceField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vector4CyclicInterfaceField::vector4CyclicInterfaceField
(
    const lduInterface& interface,
    const cyclicLduInterfaceField& field
)
:
    LduInterfaceField<vector4>(interface),
    cyclicInterface_(refCast<const cyclicLduInterface>(interface)),
    doTransform_(field.doTransform()),
    rank_(field.rank())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vector4CyclicInterfaceField::updateInterfaceMatrix
(
    Field<vector4>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<vector4>& psiInternal,
    const Field<tensor4>& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFc =
        lduAddr.patchAddr(cyclicInterface_.neighbPatchID());

    Field<vector4> pnf(nbrFc.size());

    foamExecutor exec;
    {
        vector4* __restrict__ pnfPtr = pnf.begin();
        const vector4* const __restrict__ psiPtr = psiInternal.cbegin();
        const label* const __restrict__ nbrFcPtr = nbrFc.cbegin();

        auto Lambda = [=](label facei)
        {
            pnfPtr[facei] = psiPtr[nbrFcPtr[facei]];
        };
        exec.parallelFor(Lambda, nbrFc.size());
    }

    transformCoupleField(pnf);

    const labelUList& fc = lduAddr.patchAddr(patchId);

    vector4* __restrict__ resultPtr = result.begin();
    const tensor4* const __restrict__ coeffsPtr = coeffs.cbegin();
    const vector4* const __restrict__ pnfPtr = pnf.cbegin();
    const label* const __restrict__ fcPtr = fc.cbegin();

    // Interface coefficients are stored as actual matrix entries
    if (add)
    {
        auto Lambda = [=](label facei)
        {
            foamAtomic::AtomicAdd
            (
                resultPtr[fcPtr[facei]],
                dot(coeffsPtr[facei], pnfPtr[facei])
            );
        };
        exec.parallelFor(Lambda, fc.size());
    }
    else
    {
        auto Lambda = [=](label facei)
        {
            foamAtomic::AtomicAdd
            (
                resultPtr[fcPtr[facei]],
                -dot(coeffsPtr[facei], pnfPtr[facei])
            );
        };
        exec.parallelFor(Lambda, fc.size());
    }
}


// ************************************************************************* //
