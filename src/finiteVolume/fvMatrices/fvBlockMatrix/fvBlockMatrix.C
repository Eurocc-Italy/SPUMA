/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2022 Stefano Oliani
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

#include "fvBlockMatrix.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class sourceType, class blockType>
template<class Type2>
void Foam::fvBlockMatrix<sourceType,blockType>::addToInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    foamExecutor exec;
    const Type2* const __restrict__ pfPtr = pf.cbegin();
    Type2* __restrict__ intfPtr = intf.begin();
    const label* const __restrict__ addrPtr = addr.cbegin();

    auto LambdaAdd = [=](label facei)
    {
        foamAtomic::AtomicAdd(intfPtr[addrPtr[facei]], pfPtr[facei]);
    };
    exec.parallelFor(LambdaAdd, addr.size());
}

template<class sourceType, class blockType>
template<class Type2>
void Foam::fvBlockMatrix<sourceType,blockType>::addToInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    addToInternalField(addr, tpf(), intf);
    tpf.clear();
}

template<class sourceType, class blockType>
template<class Type2>
void Foam::fvBlockMatrix<sourceType,blockType>::subtractFromInternalField
(
    const labelUList& addr,
    const Field<Type2>& pf,
    Field<Type2>& intf
) const
{
    if (addr.size() != pf.size())
    {
        FatalErrorInFunction
            << "sizes of addressing and field are different"
            << abort(FatalError);
    }

    foamExecutor exec;

    const Type2* const __restrict__ pfPtr = pf.cbegin();
    Type2* __restrict__ intfPtr = intf.begin();
    const label* const __restrict__ addrPtr = addr.cbegin();

    auto LambdaAdd = [=](label facei)
    {
        foamAtomic::AtomicAdd(intfPtr[addrPtr[facei]], -pfPtr[facei]);
    };
    exec.parallelFor(LambdaAdd, addr.size());
}

template<class sourceType, class blockType>
template<class Type2>
void Foam::fvBlockMatrix<sourceType,blockType>::subtractFromInternalField
(
    const labelUList& addr,
    const tmp<Field<Type2>>& tpf,
    Field<Type2>& intf
) const
{
    subtractFromInternalField(addr, tpf(), intf);
    tpf.clear();
}

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType,blockType>::addBoundaryDiag
(
    fvMatrix<sourceType> matrix,
    scalarField& diag,
    const direction solveCmpt
) const
{
    forAll(matrix.internalCoeffs(), patchi)
    {
        addToInternalField
        (
            this->lduAddr().patchAddr(patchi),
            matrix.internalCoeffs()[patchi].component(solveCmpt),
            diag
        );
    }
}

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType,blockType>::addBoundarySource
(
    fvMatrix<sourceType> matrix,
    Field<sourceType>& source,
    const bool couples
) const
{
    forAll(matrix.psi().boundaryField(), patchi)
    {
        const fvPatchField<sourceType>& ptf = matrix.psi().boundaryField()[patchi];
        const Field<sourceType>& pbc = matrix.boundaryCoeffs()[patchi];

        if (!ptf.coupled())
        {
            addToInternalField(this->lduAddr().patchAddr(patchi), pbc, source);
        }
        else if (couples)
        {
            const tmp<Field<sourceType>> tpnf = ptf.patchNeighbourField();
            const Field<sourceType>& pnf = tpnf();

            const labelUList& addr = this->lduAddr().patchAddr(patchi);

            foamExecutor exec;

            const label* const __restrict__ addrPtr = addr.cbegin();
            const sourceType* const __restrict__ pbcPtr = pbc.cbegin();
            const sourceType* const __restrict__ pnfPtr = pnf.cbegin();
            sourceType* __restrict__ sourcePtr = source.begin();

            auto LambdaAdd = [=](label facei)
            {
                foamAtomic::AtomicAdd
                (
                    sourcePtr[addrPtr[facei]],
                    cmptMultiply(pbcPtr[facei], pnfPtr[facei])
                );
            };
            exec.parallelFor(LambdaAdd, addr.size());
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType,blockType>::fvBlockMatrix
(
    const lduMesh& mesh
)
:
    refCount(),
    LduMatrix<sourceType,blockType,blockType>(mesh)
{}

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType,blockType>::fvBlockMatrix
(
    const fvBlockMatrix<sourceType,blockType>& bfm
)
:
    refCount(),
    LduMatrix<sourceType,blockType,blockType>(bfm)
{
    // Perform assignment of interfacesUpper and interfacesLower, which doesn't happen in LduMatrix
    this->interfacesUpper().resize(bfm.interfacesUpper().size());
    this->interfacesLower().resize(bfm.interfacesUpper().size());

    forAll(this->interfacesUpper(), i)
    {
        this->interfacesUpper().set(i, new Field<blockType>(bfm.interfacesUpper()[i]));
    }
    forAll(this->interfacesLower(), i)
    {
        this->interfacesLower().set(i, new Field<blockType>(bfm.interfacesLower()[i]));
    }
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType,blockType>> Foam::fvBlockMatrix<sourceType,blockType>::clone() const
{
    return tmp<fvBlockMatrix<sourceType,blockType>>
    (
        new fvBlockMatrix<sourceType,blockType>(*this)
    );
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType,blockType>::~fvBlockMatrix()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType,blockType>::insertBlock
(
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& leftField,
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& rightField
)
{
    const fvMesh& mesh = leftField.mesh();
    const labelUList& owner = mesh.owner();

    fvBlockMatrix<sourceType, blockType> mx(mesh);

    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();
    Field<blockType>& diag = mx.diag();

    foamExecutor exec;

    blockType* __restrict__ uppPtr = upp.begin();
    blockType* __restrict__ lowPtr = low.begin();
    blockType* __restrict__ diagPtr = diag.begin();
    const blockType* const __restrict__ rightFieldPtr = rightField.primitiveField().cbegin();
    const blockType* const __restrict__ leftFieldPtr = leftField.primitiveField().cbegin();
    const scalar* const __restrict__ magSfPtr = mesh.magSf().cbegin();

    auto LambdaInsertBlock = [=](label facei)
    {
        const scalar magSfi = magSfPtr[facei];
        uppPtr[facei] = 0.5 * magSfi * rightFieldPtr[facei];
        lowPtr[facei] = -0.5 * magSfi * leftFieldPtr[facei];
    };
    exec.parallelFor(LambdaInsertBlock, owner.size());

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];

        const Field<blockType>& leftFp = leftField.boundaryField()[patchi];
        const Field<blockType>& rightFp = rightField.boundaryField()[patchi];

        mx.interfacesUpper().set(patchi, 0.5*magSfp * rightFp);
        mx.interfacesLower().set(patchi, -0.5*magSfp * leftFp);

        // Don't include physical boundaries because those are dealt with separately
        if (mesh.boundary()[patchi].coupled())
        {
            const labelUList& bOwner = mesh.boundary()[patchi].faceCells();
            Field<blockType>& low = mx.interfacesLower()[patchi];

            const label* const __restrict__ bOwnerPtr = bOwner.cbegin();
            const blockType* const __restrict__ lowPtr = low.cbegin();

            auto LambdaInsertPatchBlock = [=](label facei)
            {
                foamAtomic::AtomicAdd(diagPtr[bOwnerPtr[facei]], -lowPtr[facei]);
            };
            exec.parallelFor(LambdaInsertPatchBlock, bOwner.size());
        }
    }

    operator+=(mx);
}

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType,blockType>::insertDissipationBlock
(
    const GeometricField<blockType, fvsPatchField, surfaceMesh>& dissField
)
{
    const fvMesh& mesh = dissField.mesh();
    const labelUList& owner = mesh.owner();

    fvBlockMatrix<sourceType, blockType> mx(mesh);

    Field<blockType>& upp = mx.upper();
    Field<blockType>& low = mx.lower();
    Field<blockType>& diag = mx.diag();

    forAll(owner, facei)
    {
        scalar magSfi = mesh.magSf()[facei];

        upp[facei] = 0.5 * magSfi * dissField[facei];
        low[facei] = 0.5 * magSfi * dissField[facei];
    }

    mx.negSumDiag();

    mx.interfacesUpper().resize(mesh.boundary().size());
    mx.interfacesLower().resize(mesh.boundary().size());

    forAll(mesh.boundary(), patchi)
    {
        const scalarField& magSfp = mesh.magSf().boundaryField()[patchi];
        const Field<blockType>& dissFp = dissField.boundaryField()[patchi];

        mx.interfacesUpper().set(patchi, 0.5*magSfp * dissFp);
        mx.interfacesLower().set(patchi, 0.5*magSfp * dissFp);

        // Include physical boundaries or not?
        // if (mesh.boundary()[patchi].coupled())
        // {
            const labelUList& bOwner = mesh.boundary()[patchi].faceCells();

            Field<blockType>& low = mx.interfacesLower()[patchi];

            forAll(bOwner, facei)
            {
                const label own = bOwner[facei];

                diag[own] -= low[facei];
            }
        // }
    }

    operator-=(mx);
}


template<class sourceType, class blockType>
template<class psiType>
void Foam::fvBlockMatrix<sourceType,blockType>::Amul
(
    Field<sourceType>& Apsi,
    const GeometricField<psiType, fvPatchField, volMesh>& psi,
    const fvMesh& mesh
) const
{
    scalarField interfaceDummyCoeffs(1);

    sourceType* __restrict__ ApsiPtr = Apsi.begin();
    const psiType* const __restrict__ psiPtr = psi.primitiveField().begin();
    const blockType* const __restrict__ diagPtr = this->diag().begin();

    const label* const __restrict__ uPtr = this->lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = this->lduAddr().lowerAddr().begin();

    LduInterfaceFieldPtrsList<psiType> interfaces =
    psi.boundaryField().interfaces();

    Field<psiType> dummyResult(1);

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initProcessorInterfaces
    (
        psi.primitiveField(),
        dummyResult,
        interfaces
    );

    foamExecutor exec;

    const label nCells = this->diag().size();
    auto LamdaDiag = [=](label cell)
    {
        ApsiPtr[cell] = dot(diagPtr[cell], psiPtr[cell]);
    };
    exec.parallelFor(LamdaDiag, nCells);

    if (this->hasLower() || this->hasUpper())
    {
        const blockType* const __restrict__ upperPtr = this->upper().begin();
        const blockType* const __restrict__ lowerPtr = this->lower().begin();

        const label nFaces = this->upper().size();

        auto LambdaOffDiag = [=](label face)
        {
            foamAtomic::AtomicAdd(ApsiPtr[uPtr[face]], dot(lowerPtr[face], psiPtr[lPtr[face]]));
            foamAtomic::AtomicAdd(ApsiPtr[lPtr[face]], dot(upperPtr[face], psiPtr[uPtr[face]]));
        };
        exec.parallelFor(LambdaOffDiag, nFaces);
    }

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
    )
    {
        forAll(mesh.boundary(), patchi)
        {
            if (isA<processorFvPatch>(mesh.boundary()[patchi]))
            {
                const processorFvPatchField<psiType>& procField =
                    (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[patchi]));

                updateProcessorInterfaces
                (
                    Apsi,
                    Pstream::defaultCommsType,
                    procField,
                    mesh,
                    patchi
                );
            }
            else if (mesh.boundary()[patchi].coupled())
            {
                const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

                // Diagonal matrix does not have interfaces defined
                if (this->interfacesUpper().size())
                {
                    Field<psiType> pdS = psi.boundaryField()[patchi].patchNeighbourField();
                    const Field<blockType>& pdSByS = this->interfacesUpper()[patchi];

                    const blockType* const __restrict__ pdSBySPtr = pdSByS.cbegin();
                    const psiType* const __restrict__ pdSPtr = pdS.cbegin();
                    const label* const __restrict__ faceCellsPtr = faceCells.cbegin();

                    auto LambdaOffProc = [=](label facei)
                    {
                        foamAtomic::AtomicAdd(ApsiPtr[faceCellsPtr[facei]], dot(pdSBySPtr[facei], pdSPtr[facei]));
                    };
                    exec.parallelFor(LambdaOffProc, pdS.size());
                }

                if(psi.boundaryField()[patchi].coupledDeleteRecv())
                {
                    // if (debug)
                    // {
                    //     Info << "deleting Recv requests (ugly)"
                    // }
                }
            }
        }
    }
    else if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        // Try and consume interfaces as they become available
        bool allUpdated = false;

        // Block for everything
        if (Pstream::parRun())
        {
            if (allUpdated)
            {
                // All received. Just remove all outstanding requests
                UPstream::resetRequests(startRequest);
            }
            else
            {
                // Block for all requests and remove storage
                UPstream::waitRequests(startRequest);
            }
        }

        // Consume
        forAll(interfaces, interfacei)
        {
            if
            (
                interfaces.set(interfacei)
            && !interfaces[interfacei].updatedMatrix()
            )
            {
                if (isA<processorFvPatch>(mesh.boundary()[interfacei]))
                {
                    const processorFvPatchField<psiType>& procField =
                        (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

                    updateProcessorInterfaces
                    (
                        Apsi,
                        Pstream::defaultCommsType,
                        procField,
                        mesh,
                        interfacei
                    );
                }
                else if (mesh.boundary()[interfacei].coupled())
                {
                    const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

                    // Diagonal matrix does not have interfaces defined
                    if (this->interfacesUpper().size())
                    {
                        Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();
                        const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

                        const blockType* const __restrict__ pdSBySPtr = pdSByS.cbegin();
                        const psiType* const __restrict__ pdSPtr = pdS.cbegin();
                        const label* const __restrict__ faceCellsPtr = faceCells.cbegin();

                        auto LambdaOffProc = [=](label facei)
                        {
                            foamAtomic::AtomicAdd(ApsiPtr[faceCellsPtr[facei]], dot(pdSBySPtr[facei], pdSPtr[facei]));
                        };
                        exec.parallelFor(LambdaOffProc, pdS.size());
                    }
                    if(psi.boundaryField()[interfacei].coupledDeleteRecv())
                    {
                        // if (debug)
                        // {
                        //     Info << "deleting Recv requests (ugly)"
                        // }
                    }
                }
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll(patchSchedule, i)
        {
            label interfacei = patchSchedule[i].patch;

            if (mesh.boundary()[interfacei].coupled())
            {
                if (patchSchedule[i].init)
                {
                    interfaces[interfacei].initInterfaceMatrixUpdate
                    (
                        dummyResult,
                        true,
                        this->lduAddr(),
                        interfacei,
                        psi.primitiveField(),
                        interfaceDummyCoeffs,
                        Pstream::defaultCommsType
                    );
                }
                else
                {
                    if(isA<processorFvPatch>(mesh.boundary()[interfacei]))
                    {
                        const processorFvPatchField<psiType>& procField =
                            (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

                        updateProcessorInterfaces
                        (
                            Apsi,
                            Pstream::defaultCommsType,
                            procField,
                            mesh,
                           interfacei
                        );
                    }
                    else
                    {
                        const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

                        // Diagonal matrix does not have interfaces defined
                        if (this->interfacesUpper().size())
                        {
                            Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();
                            const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

                            const blockType* const __restrict__ pdSBySPtr = pdSByS.cbegin();
                            const psiType* const __restrict__ pdSPtr = pdS.cbegin();
                            const label* const __restrict__ faceCellsPtr = faceCells.cbegin();

                            auto LambdaOffProc = [=](label facei)
                            {
                                foamAtomic::AtomicAdd
                                (
                                    ApsiPtr[faceCellsPtr[facei]], dot(pdSBySPtr[facei], pdSPtr[facei])
                                );
                            };
                            exec.parallelFor(LambdaOffProc, pdS.size());
                        }
                        if(psi.boundaryField()[interfacei].coupledDeleteRecv())
                        {
                            // if (debug)
                            // {
                            //     Info << "deleting Recv requests (ugly)"
                            // }
                        }
                    }
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfacei=patchSchedule.size()/2;
            interfacei<interfaces.size();
            interfacei++
        )
        {
            if(isA<processorFvPatch>(mesh.boundary()[interfacei]))
            {
               const processorFvPatchField<psiType>& procField =
                   (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

               updateProcessorInterfaces
               (
                   Apsi,
                   Pstream::defaultCommsType,
                   procField,
                   mesh,
                   interfacei
               );
            }
            else if (mesh.boundary()[interfacei].coupled())
            {
               const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

               // Diagonal matrix does not have interfaces defined
               if (this->interfacesUpper().size())
               {
                   Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();
                   const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

                   const blockType* const __restrict__ pdSBySPtr = pdSByS.cbegin();
                   const psiType* const __restrict__ pdSPtr = pdS.cbegin();
                   const label* const __restrict__ faceCellsPtr = faceCells.cbegin();

                   auto LambdaOffProc = [=](label facei)
                   {
                       foamAtomic::AtomicAdd(ApsiPtr[faceCellsPtr[facei]], dot(pdSBySPtr[facei], pdSPtr[facei]));
                   };
                   exec.parallelFor(LambdaOffProc, pdS.size());
                }
                if(psi.boundaryField()[interfacei].coupledDeleteRecv())
                {
                    // if (debug)
                    // {
                    //     Info << "deleting Recv requests (ugly)"
                    // }
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}

template<class sourceType, class blockType>
template<class psiType>
void Foam::fvBlockMatrix<sourceType,blockType>::AmulNoDiag
(
    Field<sourceType>& Apsi,
    const GeometricField<psiType, fvPatchField, volMesh>& psi,
    const fvMesh& mesh
) const
{
    Apsi = Zero;
    scalarField interfaceDummyCoeffs(1);

    sourceType* __restrict__ ApsiPtr = Apsi.begin();

    const psiType* const __restrict__ psiPtr = psi.primitiveField().begin();

    const label* const __restrict__ uPtr = this->lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = this->lduAddr().lowerAddr().begin();

    LduInterfaceFieldPtrsList<psiType> interfaces =
        psi.boundaryField().interfaces();

    Field<psiType> dummyResult(1);

    const label startRequest = Pstream::nRequests();

    // Initialise the update of interfaced interfaces
    initProcessorInterfaces
    (
        psi.primitiveField(),
        dummyResult,
        interfaces
    );

    if (this->hasLower() || this->hasUpper())
    {
        const blockType* const __restrict__ upperPtr = this->upper().begin();
        const blockType* const __restrict__ lowerPtr = this->lower().begin();

        const label nFaces = this->upper().size();

        for (label face=0; face<nFaces; face++)
        {
            ApsiPtr[uPtr[face]] += dot(lowerPtr[face], psiPtr[lPtr[face]]);
            ApsiPtr[lPtr[face]] += dot(upperPtr[face], psiPtr[uPtr[face]]);
        }
    }

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
    )
    {
        forAll(mesh.boundary(), patchi)
        {
            if (isA<processorFvPatch>(mesh.boundary()[patchi]))
            {
                const processorFvPatchField<psiType>& procField =
                    (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[patchi]));

                updateProcessorInterfaces
                (
                    Apsi,
                    Pstream::defaultCommsType,
                    procField,
                    mesh,
                    patchi
                );
            }
            else if (mesh.boundary()[patchi].coupled())
            {
                const labelUList& faceCells = mesh.boundary()[patchi].faceCells();

                // Diagonal matrix does not have interfaces defined
                if (this->interfacesUpper().size())
                {
                    Field<psiType> pdS = psi.boundaryField()[patchi].patchNeighbourField();

                    const Field<blockType>& pdSByS = this->interfacesUpper()[patchi];

                    forAll(pdS, facei)
                    {
                        Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
                    }
                }
                if(psi.boundaryField()[patchi].coupledDeleteRecv())
                {
                    // if (debug)
                    // {
                    //     Info << "deleting Recv requests (ugly)"
                    // }
                }
            }
        }
    }
    else if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        // Try and consume interfaces as they become available
        bool allUpdated = false;

        // Block for everything
        if (Pstream::parRun())
        {
            if (allUpdated)
            {
                // All received. Just remove all outstanding requests
                UPstream::resetRequests(startRequest);
            }
            else
            {
                // Block for all requests and remove storage
                UPstream::waitRequests(startRequest);
            }
        }

        // Consume
        forAll(interfaces, interfacei)
        {
            if
            (
                interfaces.set(interfacei)
            && !interfaces[interfacei].updatedMatrix()
            )
            {
                if (isA<processorFvPatch>(mesh.boundary()[interfacei]))
                {
                    const processorFvPatchField<psiType>& procField =
                        (refCast<const processorFvPatchField<psiType>>(psi.boundaryField()[interfacei]));

                    updateProcessorInterfaces
                    (
                        Apsi,
                        Pstream::defaultCommsType,
                        procField,
                        mesh,
                        interfacei
                    );
                }
                else if (mesh.boundary()[interfacei].coupled())
                {
                    const labelUList& faceCells = mesh.boundary()[interfacei].faceCells();

                    // Diagonal matrix does not have interfaces defined
                    if (this->interfacesUpper().size())
                    {
                        Field<psiType> pdS = psi.boundaryField()[interfacei].patchNeighbourField();

                        const Field<blockType>& pdSByS = this->interfacesUpper()[interfacei];

                        forAll(pdS, facei)
                        {
                            Apsi[faceCells[facei]] += dot(pdSByS[facei],pdS[facei]);
                        }
                    }
                    if(psi.boundaryField()[interfacei].coupledDeleteRecv())
                    {
                        // if (debug)
                        // {
                        //     Info << "deleting Recv requests (ugly)"
                        // }
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType,blockType>::insertEquation
(
    fvMatrix<sourceType>& matrix
)
{}

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType, blockType>::operator=
(
    const fvBlockMatrix<sourceType, blockType>& blfm
)
{
    LduMatrix<sourceType, blockType, blockType>::operator=(std::move(blfm));

    // Perform assignment of interfacesUpper and interfacesLower, which doesn't happen in LduMatrix
    this->interfacesUpper().resize(blfm.interfacesUpper().size());
    this->interfacesLower().resize(blfm.interfacesUpper().size());

    forAll(this->interfacesUpper(), i)
    {
        this->interfacesUpper().set(i, new Field<blockType>(blfm.interfacesUpper()[i]));
    }
    forAll(this->interfacesLower(), i)
    {
        this->interfacesLower().set(i, new Field<blockType>(blfm.interfacesLower()[i]));
    }
}

template<class sourceType, class blockType>
void Foam::fvBlockMatrix<sourceType, blockType>::operator=
(
    const tmp<fvBlockMatrix<sourceType, blockType> >& tblfm
)
{
    operator=(tblfm());
    tblfm.clear();
}

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType, blockType>& Foam::fvBlockMatrix<sourceType, blockType>::operator+=
(
    const fvBlockMatrix<sourceType,
    blockType>& B
)
{
    // If B has extra interfacesUpper or interfacesLower not present here,
    // we have to expand to accommodate them before calling LduMatrix::operator+=
    const label nBupper = B.interfacesUpper().size();
    const label nupper = this->interfacesUpper().size();
    if (nBupper > nupper)
    {
        this->interfacesUpper().resize(nBupper);
        for(label i = nupper; i < nBupper; i++)
        {
            this->interfacesUpper().set(i,new Field<blockType>(B.interfacesUpper()[i].size(), pTraits<blockType>::zero));
        }
    }
    const label nBlower = B.interfacesLower().size();
    const label nlower = this->interfacesLower().size();
    if (nBlower > nlower)
    {
        this->interfacesLower().resize(nBlower);
        for(label i = nlower; i < nBlower; i++)
        {
            this->interfacesLower().set(i, new Field<blockType>(B.interfacesLower()[i].size(), pTraits<blockType>::zero));
        }
    }

    LduMatrix<sourceType, blockType, blockType>::operator+=(B);

    return *this;
}

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType, blockType>& Foam::fvBlockMatrix<sourceType, blockType>::operator+=
(
    const tmp<fvBlockMatrix<sourceType, blockType> >& tB
)
{
    operator+=(tB());
    tB.clear();
    return *this;
}

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType, blockType>& Foam::fvBlockMatrix<sourceType, blockType>::operator-=
(
    const fvBlockMatrix<sourceType, blockType>& B
)
{
    // If B has extra interfacesUpper or interfacesLower not present here,
    // we have to expand to accommodate them before calling LduMatrix::operator+=
    const label nBupper = B.interfacesUpper().size();
    const label nupper = this->interfacesUpper().size();
    if (nBupper > nupper)
    {
        this->interfacesUpper().resize(nBupper);
        for(label i = nupper; i < nBupper; i++)
        {
            this->interfacesUpper().set
            (
                i,
                new Field<blockType>
                (
                    B.interfacesUpper()[i].size(),
                    pTraits<blockType>::zero
                )
            );
        }
    }
    const label nBlower = B.interfacesLower().size();
    const label nlower = this->interfacesLower().size();
    if (nBlower > nlower)
    {
        this->interfacesLower().resize(nBlower);
        for(label i = nlower; i < nBlower; i++)
        {
            this->interfacesLower().set
            (
                i,
                new Field<blockType>
                (
                    B.interfacesLower()[i].size(),
                    pTraits<blockType>::zero
                )
            );
        }
    }

    LduMatrix<sourceType, blockType, blockType>::operator-=(B);
    return *this;
}

template<class sourceType, class blockType>
Foam::fvBlockMatrix<sourceType, blockType>& Foam::fvBlockMatrix<sourceType, blockType>::operator-=
(
    const tmp<fvBlockMatrix<sourceType,blockType> >& tB
)
{
    operator-=(tB());
    tB.clear();
    return *this;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator+
(
    const fvBlockMatrix<sourceType, blockType>& A,
    const fvBlockMatrix<sourceType, blockType>& B
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC
    (
        new fvBlockMatrix<sourceType, blockType>(A)
    );
    tC() += B;
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator+
(
    const tmp<fvBlockMatrix<sourceType, blockType> >& tA,
    const fvBlockMatrix<sourceType, blockType>& B
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC(tA.ptr());
    tC() += B;
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator+
(
    const fvBlockMatrix<sourceType, blockType>& A,
    const tmp<fvBlockMatrix<sourceType, blockType> >& tB
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC(tB.ptr());
    tC() += A;
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator+
(
    const tmp<fvBlockMatrix<sourceType, blockType> >& tA,
    const tmp<fvBlockMatrix<sourceType, blockType> >& tB
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC(tA.ptr());
    tC.ref() += tB();
    tB.clear();
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator-
(
    const fvBlockMatrix<sourceType, blockType>& A,
    const fvBlockMatrix<sourceType, blockType>& B
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC
    (
        new fvBlockMatrix<sourceType, blockType>(A)
    );
    tC() -= B;
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator-
(
    const tmp<fvBlockMatrix<sourceType, blockType> >& tA,
    const fvBlockMatrix<sourceType, blockType>& B
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC(tA.ptr());
    tC() -= B;
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator-
(
    const fvBlockMatrix<sourceType, blockType>& A,
    const tmp<fvBlockMatrix<sourceType, blockType> >& tB
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC(tB.ptr());
    tC() -= A;
    return tC;
}

template<class sourceType, class blockType>
Foam::tmp<Foam::fvBlockMatrix<sourceType, blockType> > Foam::operator-
(
    const tmp<fvBlockMatrix<sourceType, blockType> >& tA,
    const tmp<fvBlockMatrix<sourceType, blockType> >& tB
)
{
    tmp<fvBlockMatrix<sourceType, blockType> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}

//- Multiplication by a uniform value (a primitive)

template<class Type2>
Foam::tmp<Foam::fvBlockMatrix<Foam::vector, Type2> > Foam::operator*
(
    const tmp<fvBlockMatrix<scalar, scalar> >& tA, const Type2& B
)
{
    tmp<fvBlockMatrix<vector, Type2> > tC(new fvBlockMatrix<vector, Type2>(tA->mesh()));

    if (tA->hasDiag())
    {
        tC->diag() = tA->diag()*B;
    }

    if (tA->hasUpper())
    {
        tC->upper() = tA->upper()*B;
    }

    if (tA->hasLower())
    {
        tC->lower() = tA->lower()*B;
    }

    tC->interfacesUpper() = tA->interfacesUpper()*B;
    tC->interfacesLower() = tA->interfacesLower()*B;

    tA.clear();
    return tC;
}

// ************************************************************************* //
