/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |                                                 |
     \\/     M anipulation  |                                                 |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Cineca
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_lduCalculatedCyclicField_H
#define Foam_lduCalculatedCyclicField_H

#include "LduInterfaceField.H"
#include "cyclicLduInterface.H"
#include "cyclicLduInterfaceField.H"
#include "tensor4.H"
#include "vector4.H"

namespace Foam
{

template<class Type>
class lduCalculatedCyclicField
:
    public LduInterfaceField<Type>,
    public cyclicLduInterfaceField
{
    const cyclicLduInterface& cyclicInterface_;
    bool doTransform_;
    int rank_;

public:

    lduCalculatedCyclicField
    (
        const lduInterface& interface,
        const cyclicLduInterfaceField& field
    )
    :
        LduInterfaceField<Type>(interface),
        cyclicInterface_(refCast<const cyclicLduInterface>(interface)),
        doTransform_(field.doTransform()),
        rank_(field.rank())
    {}

    virtual void updateInterfaceMatrix
    (
        Field<Type>& result,
        const bool add,
        const lduAddressing&,
        const label,
        const Field<Type>& psiInternal,
        const scalarField& coeffs,
        const Pstream::commsTypes
    ) const
    {
        const labelUList& faceCells = this->interface().faceCells();
        const labelUList& neighbourFaceCells =
            refCast<const lduInterface>(cyclicInterface_.neighbPatch())
           .faceCells();

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
        const auto& neighbourFaceCells =
            lduAddr.patchAddr(cyclicInterface_.neighbPatchID());
        solveScalarField patchNeighbourField(psiInternal, neighbourFaceCells);

        cyclicLduInterfaceField::transformCoupleField
        (
            patchNeighbourField,
            cmpt
        );

        const auto& faceCells = lduAddr.patchAddr(patchId);
        this->addToInternalField
        (
            result,
            !add,
            faceCells,
            coeffs,
            patchNeighbourField
        );
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
            const labelUList& faceCells = this->interface().faceCells();
            const labelUList& neighbourFaceCells =
                refCast<const lduInterface>(cyclicInterface_.neighbPatch())
               .faceCells();

            Field<Type> patchNeighbourField(neighbourFaceCells.size());

            forAll(patchNeighbourField, facei)
            {
                patchNeighbourField[facei] =
                    psiInternal[neighbourFaceCells[facei]];
            }

            cyclicLduInterfaceField::transformCoupleField
            (
                patchNeighbourField
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
                    result[faceCells[facei]] +=
                        dot(coeffs[facei], patchNeighbourField[facei]);
                }
                else
                {
                    result[faceCells[facei]] -=
                        dot(coeffs[facei], patchNeighbourField[facei]);
                }
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
