/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2025 OpenCFD Ltd.
    Copyright (C) 2025 Cineca
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

#include "lduCalculatedProcessorField.H"
#include "vector4.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduInterface& interface
)
:
    LduInterfaceField<Type>(interface),
    procInterface_(refCast<const processorLduInterface>(interface)),
    doTransform_(false),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduInterface& interface,
    const processorLduInterfaceField& field
)
:
    LduInterfaceField<Type>(interface),
    procInterface_(refCast<const processorLduInterface>(interface)),
    doTransform_(field.doTransform()),
    sendRequest_(-1),
    recvRequest_(-1)
{}


template<class Type>
Foam::lduCalculatedProcessorField<Type>::lduCalculatedProcessorField
(
    const lduCalculatedProcessorField<Type>& ptf
)
:
    LduInterfaceField<Type>(refCast<const lduInterface>(ptf)),
    procInterface_(ptf.procInterface_),
    doTransform_(ptf.doTransform_),
    sendRequest_(-1),
    recvRequest_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::lduCalculatedProcessorField<Type>::all_ready() const
{
    return UPstream::finishedRequestPair(recvRequest_, sendRequest_);
}


template<class Type>
bool Foam::lduCalculatedProcessorField<Type>::ready() const
{
    const bool ok = UPstream::finishedRequest(recvRequest_);
    if (ok)
    {
        recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
    }
    return ok;
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::initFieldTransfer
(
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal
) const
{
    if (!this->all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on interface "
            << abort(FatalError);
    }

    const labelUList& fc = lduAddr.patchAddr(patchId);

    sendBuf_.resize_nocopy(fc.size());
    recvBuf_.resize_nocopy(fc.size());

    foamExecutor exec;
    auto sendBufp = sendBuf_.begin();
    const auto psiInternalp = psiInternal.cbegin();
    const auto fcp = fc.cbegin();
    auto Lambda = [=](label i)
    {
        sendBufp[i] = psiInternalp[fcp[i]];
    };
    exec.parallelFor(Lambda, fc.size());

    recvRequest_ = UPstream::nRequests();
    UIPstream::read
    (
        UPstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        recvBuf_,
        procInterface_.tag(),
        procInterface_.comm()
    );

    sendRequest_ = UPstream::nRequests();
    UOPstream::write
    (
        UPstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        sendBuf_,
        procInterface_.tag(),
        procInterface_.comm()
    );

    this->updatedMatrix(false);
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::initInterfaceMatrixUpdate
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
    initFieldTransfer(lduAddr, patchId, psiInternal);
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::initInterfaceMatrixUpdate
(
    Field<Type>&,
    const bool,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<Type>& psiInternal,
    const Field<tensor4>&,
    const Pstream::commsTypes
) const
{
    initFieldTransfer(lduAddr, patchId, psiInternal);
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const lduAddressing&,
    const label,
    const Field<Type>&,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    UPstream::waitRequest(recvRequest_);
    recvRequest_ = -1;
    if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;

    const labelUList& faceCells = this->interface().faceCells();
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


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::updateInterfaceMatrix
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

        UPstream::waitRequest(recvRequest_);
        recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;

        processorLduInterfaceField::transformCoupleField(recvBuf_);

        const labelUList& faceCells = this->interface().faceCells();
        const label nFaces = min
        (
            faceCells.size(),
            min(coeffs.size(), recvBuf_.size())
        );

        for (label facei = 0; facei < nFaces; ++facei)
        {
            if (add)
            {
                result[faceCells[facei]] += dot(coeffs[facei], recvBuf_[facei]);
            }
            else
            {
                result[faceCells[facei]] -= dot(coeffs[facei], recvBuf_[facei]);
            }
        }

        this->updatedMatrix(true);
    }
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::initInterfaceMatrixUpdate
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (!this->all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on interface "
            //<< procInterface_.name()
            << abort(FatalError);
    }

    // Bypass patchInternalField since uses fvPatch addressing
    const labelUList& fc = lduAddr.patchAddr(patchId);

    {
        scalarSendBuf_.resize_nocopy(fc.size());
        scalarRecvBuf_.resize_nocopy(fc.size());
        foamExecutor exec;
        auto sSendBufp = scalarSendBuf_.begin();
        const auto psiInternalp = psiInternal.cbegin();
        const auto fcp = fc.cbegin();
        auto Lambda = [=](label i)
        {
            sSendBufp[i] = psiInternalp[fcp[i]];
        };
        exec.parallelFor(Lambda, fc.size());


    }

    recvRequest_ = UPstream::nRequests();
    UIPstream::read
    (
        UPstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarRecvBuf_,
        procInterface_.tag(),
        procInterface_.comm()
    );

    sendRequest_ = UPstream::nRequests();
    UOPstream::write
    (
        UPstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        scalarSendBuf_,
        procInterface_.tag(),
        procInterface_.comm()
    );

    this->updatedMatrix(false);
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::addToInternalField
(
    solveScalarField& result,
    const bool add,
    const scalarField& coeffs,
    const solveScalarField& vals
) const
{
    const labelUList& faceCells = this->interface().faceCells();

    foamExecutor exec;
    auto resultp = result.begin();
    const auto faceCellsp=faceCells.cbegin();
    const auto coeffsp = coeffs.cbegin();
    const auto valsp = vals.cbegin();

    if (add)
    {
        auto Lambda = [=](label elemI)
        {
            foamAtomic::AtomicAdd(resultp[faceCellsp[elemI]], coeffsp[elemI]*valsp[elemI]);
        };
        exec.parallelFor(Lambda, faceCells.size());
    }
    else
    {
        auto Lambda = [=](label elemI)
        {
            foamAtomic::AtomicAdd(resultp[faceCellsp[elemI]], -coeffsp[elemI]*valsp[elemI]);
        };
        exec.parallelFor(Lambda, faceCells.size());
    }
}


template<class Type>
void Foam::lduCalculatedProcessorField<Type>::updateInterfaceMatrix
(
    solveScalarField& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const solveScalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    {
        // Require receive data.
        // Only update the send request state.
        UPstream::waitRequest(recvRequest_); recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
    }

    // Consume straight from receive buffer. Note use of our own
    // helper to avoid using fvPatch addressing
    addToInternalField(result, !add, coeffs, scalarRecvBuf_);

    this->updatedMatrix(true);
}


// ************************************************************************* //
