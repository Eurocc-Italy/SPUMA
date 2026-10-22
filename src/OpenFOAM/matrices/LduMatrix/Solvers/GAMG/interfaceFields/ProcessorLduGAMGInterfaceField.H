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

#ifndef Foam_ProcessorLduGAMGInterfaceField_H
#define Foam_ProcessorLduGAMGInterfaceField_H

#include "LduInterfaceField.H"
#include "processorLduGAMGInterface.H"
#include "transformField.H"
#include "processorLduInterfaceField.H"

namespace Foam
{

template<class Type>
class ProcessorLduGAMGInterfaceField
:
    public LduInterfaceField<Type>,
    public processorLduInterfaceField
{
    const processorLduGAMGInterface& procInterface_;
    bool doTransform_;
    int rank_;
    mutable label sendRequest_;
    mutable label recvRequest_;
    mutable Field<Type> sendBuf_;
    mutable Field<Type> recvBuf_;
    mutable solveScalarField scalarSendBuf_;
    mutable solveScalarField scalarRecvBuf_;

    //- Face cells in the memory pool, for the device update of block
    //- coefficients (built on first use)
    mutable labelList faceCellsPool_;

    void setPoolAddressing() const
    {
        const labelUList& faceCells = procInterface_.faceCells();

        if (faceCellsPool_.size() != faceCells.size())
        {
            faceCellsPool_ = labelList(faceCells.size(), poolSwitch(1));
            faceCellsPool_.deepCopy(faceCells);
        }
    }

public:

    ProcessorLduGAMGInterfaceField
    (
        const processorLduGAMGInterface& coarseInterface,
        const LduInterfaceField<Type>& fineInterface
    )
    :
        LduInterfaceField<Type>(coarseInterface),
        procInterface_(coarseInterface),
        doTransform_(false),
        rank_(pTraits<Type>::rank),
        sendRequest_(-1),
        recvRequest_(-1),
        faceCellsPool_(poolSwitch(1))
    {
        const auto& processorField =
            refCast<const processorLduInterfaceField>(fineInterface);

        doTransform_ = processorField.doTransform();
        rank_ = processorField.rank();
    }

    //- Construct with explicit transform settings (used by proc agglomeration)
    ProcessorLduGAMGInterfaceField
    (
        const processorLduGAMGInterface& coarseInterface,
        const bool doTransform,
        const int rank
    )
    :
        LduInterfaceField<Type>(coarseInterface),
        procInterface_(coarseInterface),
        doTransform_(doTransform),
        rank_(rank),
        sendRequest_(-1),
        recvRequest_(-1),
        faceCellsPool_(poolSwitch(1))
    {}

    virtual bool ready() const
    {
        const bool ok =
        (
            recvRequest_ < 0 || UPstream::finishedRequest(recvRequest_)
        );

        if (ok)
        {
            recvRequest_ = -1;

            if (sendRequest_ >= 0 && UPstream::finishedRequest(sendRequest_))
            {
                sendRequest_ = -1;
            }
        }

        return ok;
    }

    virtual void initInterfaceMatrixUpdate
    (
        Field<Type>&,
        const bool,
        const lduAddressing& lduAddr,
        const label patchId,
        const Field<Type>& psiInternal,
        const scalarField&,
        const Pstream::commsTypes
    ) const
    {
        const labelUList& faceCells = procInterface_.faceCells();

        sendBuf_.resize_nocopy(faceCells.size());

        forAll(faceCells, facei)
        {
            sendBuf_[facei] = psiInternal[faceCells[facei]];
        }

        recvBuf_.resize_nocopy(sendBuf_.size());

        recvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            recvBuf_.data_bytes(),
            recvBuf_.size_bytes(),
            procInterface_.tag(),
            procInterface_.comm()
        );

        sendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            sendBuf_.cdata_bytes(),
            sendBuf_.size_bytes(),
            procInterface_.tag(),
            procInterface_.comm()
        );

        this->updatedMatrix(false);
    }

    virtual void initInterfaceMatrixUpdate
    (
        Field<Type>&,
        const bool,
        const lduAddressing&,
        const label,
        const Field<Type>& psiInternal,
        const Field<tensor4>&,
        const Pstream::commsTypes
    ) const
    {
        setPoolAddressing();

        sendBuf_.resize_nocopy(faceCellsPool_.size());
        recvBuf_.resize_nocopy(sendBuf_.size());

        foamExecutor exec;
        Type* __restrict__ sendBufPtr = sendBuf_.begin();
        const Type* const __restrict__ psiPtr = psiInternal.cbegin();
        const label* const __restrict__ fcPtr = faceCellsPool_.cbegin();

        auto Lambda = [=](label facei)
        {
            sendBufPtr[facei] = psiPtr[fcPtr[facei]];
        };
        exec.parallelFor(Lambda, sendBuf_.size());

        recvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            recvBuf_.data_bytes(),
            recvBuf_.size_bytes(),
            procInterface_.tag(),
            procInterface_.comm()
        );

        sendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            sendBuf_.cdata_bytes(),
            sendBuf_.size_bytes(),
            procInterface_.tag(),
            procInterface_.comm()
        );

        this->updatedMatrix(false);
    }

    virtual void initInterfaceMatrixUpdate
    (
        solveScalarField&,
        const bool,
        const lduAddressing&,
        const label,
        const solveScalarField& psiInternal,
        const scalarField&,
        const direction,
        const Pstream::commsTypes commsType
    ) const
    {
        procInterface_.interfaceInternalField(psiInternal, scalarSendBuf_);

        if
        (
            commsType == Pstream::commsTypes::nonBlocking
         && !UPstream::floatTransfer
        )
        {
            scalarRecvBuf_.resize_nocopy(scalarSendBuf_.size());

            recvRequest_ = UPstream::nRequests();
            UIPstream::read
            (
                Pstream::commsTypes::nonBlocking,
                procInterface_.neighbProcNo(),
                scalarRecvBuf_.data_bytes(),
                scalarRecvBuf_.size_bytes(),
                procInterface_.tag(),
                procInterface_.comm()
            );

            sendRequest_ = UPstream::nRequests();
            UOPstream::write
            (
                Pstream::commsTypes::nonBlocking,
                procInterface_.neighbProcNo(),
                scalarSendBuf_.cdata_bytes(),
                scalarSendBuf_.size_bytes(),
                procInterface_.tag(),
                procInterface_.comm()
            );
        }
        else
        {
            procInterface_.compressedSend(commsType, scalarSendBuf_);
        }

        this->updatedMatrix(false);
    }

    virtual void updateInterfaceMatrix
    (
        Field<Type>& result,
        const bool add,
        const lduAddressing& lduAddr,
        const label patchId,
        const Field<Type>&,
        const scalarField& coeffs,
        const Pstream::commsTypes
    ) const
    {
        if (this->updatedMatrix())
        {
            return;
        }

        if (recvRequest_ >= 0)
        {
            UPstream::waitRequest(recvRequest_);
            recvRequest_ = -1;
        }

        if (sendRequest_ >= 0 && UPstream::finishedRequest(sendRequest_))
        {
            sendRequest_ = -1;
        }

        processorLduInterfaceField::transformCoupleField(recvBuf_);

        const labelUList& faceCells = procInterface_.faceCells();
        const label nFaces = min
        (
            faceCells.size(),
            min(coeffs.size(), recvBuf_.size())
        );

        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (add)
            {
                result[faceCells[facei]] -= coeffs[facei]*recvBuf_[facei];
            }
            else
            {
                result[faceCells[facei]] += coeffs[facei]*recvBuf_[facei];
            }
        }

        this->updatedMatrix(true);
    }

    virtual void updateInterfaceMatrix
    (
        Field<Type>& result,
        const bool add,
        const lduAddressing&,
        const label,
        const Field<Type>&,
        const Field<tensor4>& coeffs,
        const Pstream::commsTypes
    ) const
    {
        if constexpr
        (
            std::is_same<Type, vector4>::value
        )
        {
            if (this->updatedMatrix())
            {
                return;
            }

            if (recvRequest_ >= 0)
            {
                UPstream::waitRequest(recvRequest_);
                recvRequest_ = -1;
            }

            if (sendRequest_ >= 0 && UPstream::finishedRequest(sendRequest_))
            {
                sendRequest_ = -1;
            }

            processorLduInterfaceField::transformCoupleField(recvBuf_);

            setPoolAddressing();

            const label nFaces = min
            (
                faceCellsPool_.size(),
                min(coeffs.size(), recvBuf_.size())
            );

            foamExecutor exec;
            Type* __restrict__ resultPtr = result.begin();
            const tensor4* const __restrict__ coeffsPtr = coeffs.cbegin();
            const Type* const __restrict__ recvBufPtr = recvBuf_.cbegin();
            const label* const __restrict__ fcPtr = faceCellsPool_.cbegin();

            // Coefficients are actual matrix entries, as on the fine level
            const scalar sign = add ? 1 : -1;

            auto Lambda = [=](label facei)
            {
                foamAtomic::AtomicAdd
                (
                    resultPtr[fcPtr[facei]],
                    sign*dot(coeffsPtr[facei], recvBufPtr[facei])
                );
            };
            exec.parallelFor(Lambda, nFaces);

            this->updatedMatrix(true);
        }
    }

    virtual void updateInterfaceMatrix
    (
        solveScalarField& result,
        const bool add,
        const lduAddressing& lduAddr,
        const label patchId,
        const solveScalarField&,
        const scalarField& coeffs,
        const direction cmpt,
        const Pstream::commsTypes commsType
    ) const
    {
        if (this->updatedMatrix())
        {
            return;
        }

        if
        (
            commsType == Pstream::commsTypes::nonBlocking
         && !UPstream::floatTransfer
        )
        {
            if (recvRequest_ >= 0)
            {
                UPstream::waitRequest(recvRequest_);
                recvRequest_ = -1;
            }

            if (sendRequest_ >= 0 && UPstream::finishedRequest(sendRequest_))
            {
                sendRequest_ = -1;
            }
        }
        else
        {
            scalarRecvBuf_.resize_nocopy(coeffs.size());
            procInterface_.compressedReceive(commsType, scalarRecvBuf_);
        }

        processorLduInterfaceField::transformCoupleField
        (
            scalarRecvBuf_,
            cmpt
        );

        const labelUList& faceCells = procInterface_.faceCells();
        const label nFaces = min
        (
            faceCells.size(),
            min(coeffs.size(), scalarRecvBuf_.size())
        );

        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (add)
            {
                result[faceCells[facei]] -=
                    coeffs[facei]*scalarRecvBuf_[facei];
            }
            else
            {
                result[faceCells[facei]] +=
                    coeffs[facei]*scalarRecvBuf_[facei];
            }
        }

        this->updatedMatrix(true);
    }

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
        return rank_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //