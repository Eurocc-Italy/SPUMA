/*---------------------------------------------------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#ifndef Foam_CyclicLduGAMGInterfaceField_H
#define Foam_CyclicLduGAMGInterfaceField_H

#include "LduInterfaceField.H"
#include "cyclicLduGAMGInterface.H"
#include "transformField.H"
#include "cyclicLduInterfaceField.H"

namespace Foam
{

template<class Type>
class CyclicLduGAMGInterfaceField
:
    public LduInterfaceField<Type>,
    public cyclicLduInterfaceField
{
    const cyclicLduGAMGInterface& cyclicInterface_;
    bool doTransform_;
    int rank_;

    //- Face cells and neighbour face cells in the memory pool, for the
    //- device update of block coefficients (built on first use)
    mutable labelList faceCellsPool_;
    mutable labelList neighbFaceCellsPool_;

    //- Neighbour values buffer for the block-coefficient update
    mutable Field<Type> patchNeighbourField_;

    void setPoolAddressing() const
    {
        if (faceCellsPool_.size() || !cyclicInterface_.faceCells().size())
        {
            return;
        }

        const labelUList& faceCells = cyclicInterface_.faceCells();
        const labelUList& neighbFaceCells =
            refCast<const cyclicLduGAMGInterface>
            (
                cyclicInterface_.neighbPatch()
            ).faceCells();

        faceCellsPool_ = labelList(faceCells.size(), poolSwitch(1));
        faceCellsPool_.deepCopy(faceCells);

        neighbFaceCellsPool_ = labelList(neighbFaceCells.size(), poolSwitch(1));
        neighbFaceCellsPool_.deepCopy(neighbFaceCells);
    }

public:

    CyclicLduGAMGInterfaceField
    (
        const cyclicLduGAMGInterface& coarseInterface,
        const LduInterfaceField<Type>& fineInterface
    )
    :
        LduInterfaceField<Type>(coarseInterface),
        cyclicInterface_(coarseInterface),
        doTransform_(false),
        rank_(pTraits<Type>::rank),
        faceCellsPool_(poolSwitch(1)),
        neighbFaceCellsPool_(poolSwitch(1))
    {
        const auto& cyclicField = refCast<const cyclicLduInterfaceField>
        (
            fineInterface
        );

        doTransform_ = cyclicField.doTransform();
        rank_ = cyclicField.rank();
    }

    virtual void updateInterfaceMatrix
    (
        Field<Type>& result,
        const bool add,
        const lduAddressing& lduAddr,
        const label patchId,
        const Field<Type>& psiInternal,
        const scalarField& coeffs,
        const Pstream::commsTypes
    ) const
    {
        const labelUList& faceCells = cyclicInterface_.faceCells();
        const labelUList& neighbourFaceCells =
            refCast<const cyclicLduGAMGInterface>
            (
                cyclicInterface_.neighbPatch()
            ).faceCells();

        Field<Type> patchNeighbourField(neighbourFaceCells.size());

        forAll(patchNeighbourField, facei)
        {
            patchNeighbourField[facei] = psiInternal[neighbourFaceCells[facei]];
        }

        cyclicLduInterfaceField::transformCoupleField(patchNeighbourField);

        const label nFaces = min
        (
            faceCells.size(),
            min(coeffs.size(), patchNeighbourField.size())
        );

        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (add)
            {
                result[faceCells[facei]] -=
                    coeffs[facei]*patchNeighbourField[facei];
            }
            else
            {
                result[faceCells[facei]] +=
                    coeffs[facei]*patchNeighbourField[facei];
            }
        }
    }

    virtual void updateInterfaceMatrix
    (
        Field<Type>& result,
        const bool add,
        const lduAddressing&,
        const label,
        const Field<Type>& psiInternal,
        const Field<tensor4>& coeffs,
        const Pstream::commsTypes
    ) const
    {
        if constexpr
        (
            std::is_same<Type, vector4>::value
        )
        {
            setPoolAddressing();

            Field<Type>& patchNeighbourField = patchNeighbourField_;
            patchNeighbourField.resize_nocopy(neighbFaceCellsPool_.size());

            foamExecutor exec;
            {
                Type* __restrict__ pnfPtr = patchNeighbourField.begin();
                const Type* const __restrict__ psiPtr = psiInternal.cbegin();
                const label* const __restrict__ nbrFcPtr =
                    neighbFaceCellsPool_.cbegin();

                auto Lambda = [=](label facei)
                {
                    pnfPtr[facei] = psiPtr[nbrFcPtr[facei]];
                };
                exec.parallelFor(Lambda, patchNeighbourField.size());
            }

            cyclicLduInterfaceField::transformCoupleField(patchNeighbourField);

            const label nFaces = min
            (
                faceCellsPool_.size(),
                min(coeffs.size(), patchNeighbourField.size())
            );

            Type* __restrict__ resultPtr = result.begin();
            const tensor4* const __restrict__ coeffsPtr = coeffs.cbegin();
            const Type* const __restrict__ pnfPtr =
                patchNeighbourField.cbegin();
            const label* const __restrict__ fcPtr = faceCellsPool_.cbegin();

            // Coefficients are actual matrix entries, as on the fine level
            const scalar sign = add ? 1 : -1;

            auto Lambda = [=](label facei)
            {
                foamAtomic::AtomicAdd
                (
                    resultPtr[fcPtr[facei]],
                    sign*dot(coeffsPtr[facei], pnfPtr[facei])
                );
            };
            exec.parallelFor(Lambda, nFaces);
        }
    }

    virtual void updateInterfaceMatrix
    (
        solveScalarField& result,
        const bool add,
        const lduAddressing& lduAddr,
        const label patchId,
        const solveScalarField& psiInternal,
        const scalarField& coeffs,
        const direction cmpt,
        const Pstream::commsTypes
    ) const
    {
        const labelUList& faceCells = cyclicInterface_.faceCells();
        const labelUList& neighbourFaceCells =
            refCast<const cyclicLduGAMGInterface>
            (
                cyclicInterface_.neighbPatch()
            ).faceCells();

        solveScalarField patchNeighbourField(neighbourFaceCells.size());

        forAll(patchNeighbourField, facei)
        {
            patchNeighbourField[facei] = psiInternal[neighbourFaceCells[facei]];
        }

        cyclicLduInterfaceField::transformCoupleField
        (
            patchNeighbourField,
            cmpt
        );

        const label nFaces = min
        (
            faceCells.size(),
            min(coeffs.size(), patchNeighbourField.size())
        );

        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (add)
            {
                result[faceCells[facei]] -=
                    coeffs[facei]*patchNeighbourField[facei];
            }
            else
            {
                result[faceCells[facei]] +=
                    coeffs[facei]*patchNeighbourField[facei];
            }
        }
    }

    virtual bool doTransform() const
    {
        return doTransform_;
    }

    virtual const tensorField& forwardT() const
    {
        return cyclicInterface_.forwardT();
    }

    virtual const tensorField& reverseT() const
    {
        return cyclicInterface_.reverseT();
    }

    virtual int rank() const
    {
        return rank_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //