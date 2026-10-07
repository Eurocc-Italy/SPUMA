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

#include "vector4ProcessorInterfaceField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vector4ProcessorInterfaceField::vector4ProcessorInterfaceField
(
    const lduInterface& interface,
    const processorLduInterfaceField& field
)
:
    LduInterfaceField<vector4>(interface),
    procInterface_(refCast<const processorLduInterface>(interface)),
    doTransform_(field.doTransform()),
    sendRequest_(-1),
    recvRequest_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vector4ProcessorInterfaceField::all_ready() const
{
    return UPstream::finishedRequestPair(recvRequest_, sendRequest_);
}


bool Foam::vector4ProcessorInterfaceField::ready() const
{
    const bool ok = UPstream::finishedRequest(recvRequest_);
    if (ok)
    {
        recvRequest_ = -1;
        if (UPstream::finishedRequest(sendRequest_)) sendRequest_ = -1;
    }
    return ok;
}


void Foam::vector4ProcessorInterfaceField::initInterfaceMatrixUpdate
(
    Field<vector4>&,
    const bool,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<vector4>& psiInternal,
    const Field<tensor4>&,
    const Pstream::commsTypes
) const
{
    if (!all_ready())
    {
        FatalErrorInFunction
            << "Outstanding request(s) on interface " << patchId
            << abort(FatalError);
    }

    const labelUList& fc = lduAddr.patchAddr(patchId);

    sendBuf_.resize_nocopy(fc.size());
    recvBuf_.resize_nocopy(fc.size());

    foamExecutor exec;
    vector4* __restrict__ sendBufPtr = sendBuf_.begin();
    const vector4* const __restrict__ psiPtr = psiInternal.cbegin();
    const label* const __restrict__ fcPtr = fc.cbegin();

    auto Lambda = [=](label facei)
    {
        sendBufPtr[facei] = psiPtr[fcPtr[facei]];
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


void Foam::vector4ProcessorInterfaceField::updateInterfaceMatrix
(
    Field<vector4>& result,
    const bool add,
    const lduAddressing& lduAddr,
    const label patchId,
    const Field<vector4>&,
    const Field<tensor4>& coeffs,
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

    transformCoupleField(recvBuf_);

    const labelUList& fc = lduAddr.patchAddr(patchId);

    foamExecutor exec;
    vector4* __restrict__ resultPtr = result.begin();
    const tensor4* const __restrict__ coeffsPtr = coeffs.cbegin();
    const vector4* const __restrict__ recvBufPtr = recvBuf_.cbegin();
    const label* const __restrict__ fcPtr = fc.cbegin();

    // Interface coefficients are stored as actual matrix entries
    if (add)
    {
        auto Lambda = [=](label facei)
        {
            foamAtomic::AtomicAdd
            (
                resultPtr[fcPtr[facei]],
                dot(coeffsPtr[facei], recvBufPtr[facei])
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
                -dot(coeffsPtr[facei], recvBufPtr[facei])
            );
        };
        exec.parallelFor(Lambda, fc.size());
    }

    this->updatedMatrix(true);
}


// ************************************************************************* //
