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
    Foam::vector4ProcessorInterfaceField

Description
    Processor interface of the coupled vector4/tensor4 block system:
    exchanges vector4 values and applies the tensor4 interface coefficients.

SourceFiles
    vector4ProcessorInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vector4ProcessorInterfaceField_H
#define Foam_vector4ProcessorInterfaceField_H

#include "LduInterfaceField.H"
#include "processorLduInterface.H"
#include "processorLduInterfaceField.H"
#include "vector4.H"
#include "tensor4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class vector4ProcessorInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class vector4ProcessorInterfaceField
:
    public LduInterfaceField<vector4>,
    public processorLduInterfaceField
{
    // Private Data

        const processorLduInterface& procInterface_;

        bool doTransform_;

        mutable label sendRequest_;

        mutable label recvRequest_;

        mutable Field<vector4> sendBuf_;

        mutable Field<vector4> recvBuf_;


    // Private Member Functions

        bool all_ready() const;


public:

    // Constructors

        //- Construct from interface, taking the transform from field
        vector4ProcessorInterfaceField
        (
            const lduInterface& interface,
            const processorLduInterfaceField& field
        );


    //- Destructor
    virtual ~vector4ProcessorInterfaceField() = default;


    // Member Functions

        virtual bool ready() const;

        using LduInterfaceField<vector4>::initInterfaceMatrixUpdate;
        using LduInterfaceField<vector4>::updateInterfaceMatrix;

        virtual void initInterfaceMatrixUpdate
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


    // Processor coupled interface functions

        virtual label comm() const
        {
            return procInterface_.comm();
        }

        virtual int myProcNo() const
        {
            return procInterface_.myProcNo();
        }

        virtual int neighbProcNo() const
        {
            return procInterface_.neighbProcNo();
        }

        virtual bool doTransform() const
        {
            return doTransform_;
        }

        virtual const tensorField& forwardT() const
        {
            return procInterface_.forwardT();
        }

        virtual int rank() const
        {
            return pTraits<vector4>::rank;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
