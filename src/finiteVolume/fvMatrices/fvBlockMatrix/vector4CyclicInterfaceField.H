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

Class
    Foam::vector4CyclicInterfaceField

Description
    Cyclic interface of the coupled vector4/tensor4 block system:
    applies the tensor4 interface coefficients to the neighbour values.

SourceFiles
    vector4CyclicInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vector4CyclicInterfaceField_H
#define Foam_vector4CyclicInterfaceField_H

#include "LduInterfaceField.H"
#include "cyclicLduInterface.H"
#include "cyclicLduInterfaceField.H"
#include "vector4.H"
#include "tensor4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class vector4CyclicInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class vector4CyclicInterfaceField
:
    public LduInterfaceField<vector4>,
    public cyclicLduInterfaceField
{
    // Private Data

        const cyclicLduInterface& cyclicInterface_;

        bool doTransform_;

        int rank_;


public:

    // Constructors

        //- Construct from interface, taking the transform from field
        vector4CyclicInterfaceField
        (
            const lduInterface& interface,
            const cyclicLduInterfaceField& field
        );


    //- Destructor
    virtual ~vector4CyclicInterfaceField() = default;


    // Member Functions

        using LduInterfaceField<vector4>::updateInterfaceMatrix;

        virtual void updateInterfaceMatrix
        (
            Field<vector4>& result,
            const bool add,
            const lduAddressing& lduAddr,
            const label patchId,
            const Field<vector4>& psiInternal,
            const Field<tensor4>& coeffs,
            const Pstream::commsTypes commsType
        ) const;

        virtual void updateInterfaceMatrix
        (
            Field<vector4>&,
            const bool,
            const lduAddressing&,
            const label,
            const Field<vector4>&,
            const scalarField&,
            const Pstream::commsTypes
        ) const
        {
            NotImplemented;
        }

        virtual void updateInterfaceMatrix
        (
            solveScalarField&,
            const bool,
            const lduAddressing&,
            const label,
            const solveScalarField&,
            const scalarField&,
            const direction,
            const Pstream::commsTypes
        ) const
        {
            NotImplemented;
        }


    // Cyclic coupled interface functions

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
