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

\\*---------------------------------------------------------------------------*/

#include "LduGAMGSolver.H"
#include "LduGAMGInterface.H"
#include "LduGAMGInterfaceField.H"
#include "processorLduGAMGInterface.H"
#include "ProcessorLduGAMGInterfaceField.H"
#include "PstreamBuffers.H"
#include "OCharStream.H"
#include "ISpanStream.H"
#include "bitSet.H"
#include "SubField.H"
#include "SubList.H"
#include "scalarMatrices.H"
#include "processorLduInterface.H"
#include "cyclicLduInterface.H"
#include "tensor4.H"
#include "vector4.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduGAMGSolver<Type, DType, LUType>::LduGAMGSolver
(
    const word& fieldName,
    const LduMatrix<Type, DType, LUType>& matrix,
    const dictionary& solverControls
)
:
    LduMatrix<Type, DType, LUType>::solver
    (
        fieldName,
        matrix,
        solverControls
    ),

    // Default values for all controls
    // which may be overridden by those in controlDict
    nPreSweeps_(0),
    preSweepsLevelMultiplier_(1),
    maxPreSweeps_(4),
    nPostSweeps_(2),
    postSweepsLevelMultiplier_(1),
    maxPostSweeps_(4),
    nFinestSweeps_(2),

    interpolateCorrection_(false),
    scaleCorrection_(matrix.symmetric()),
    directSolveCoarsest_(false),

    agglomeration_
    (
        LduGAMGAgglomeration::New
        (
            matrix,
            solverControls
        )
    ),

    matrixLevels_(agglomeration_.size()),
    primitiveInterfaceLevels_(agglomeration_.size()),
    interfaceLevels_(agglomeration_.size()),
    interfaceLevelsUpper_(agglomeration_.size()),
    interfaceLevelsLower_(agglomeration_.size())
{
    readControls();

    if (agglomeration_.doProcessorAgglomerate())
    {
        forAll(agglomeration_, fineLevelIndex)
        {
            if (!agglomeration_.hasMeshLevel(fineLevelIndex))
            {
                // No mesh. Not involved in the calculation anymore
                continue;
            }

            if
            (
                (fineLevelIndex + 1) < agglomeration_.size()
             && agglomeration_.hasProcMesh(fineLevelIndex + 1)
            )
            {
                // Construct matrix without referencing the (combined) coarse
                // mesh: build a dummy interface set instead.  These coarse
                // interfaces only have to get the local matrix through
                // agglomerateMatrix; the result is overwritten by the call
                // to procAgglomerateMatrix.
                const lduInterfacePtrsList& fineMeshInterfaces =
                    agglomeration_.interfaceLevel(fineLevelIndex);

                PtrList<LduGAMGInterface> dummyPrimMeshInterfaces
                (
                    fineMeshInterfaces.size()
                );
                lduInterfacePtrsList dummyMeshInterfaces
                (
                    dummyPrimMeshInterfaces.size()
                );

                OCharStream os(IOstreamOption::BINARY);
                ISpanStream is(IOstreamOption::BINARY);

                forAll(fineMeshInterfaces, intI)
                {
                    if (fineMeshInterfaces.set(intI))
                    {
                        os.rewind();

                        refCast<const LduGAMGInterface>
                        (
                            fineMeshInterfaces[intI]
                        ).write(os);

                        is.reset(os.view());

                        dummyPrimMeshInterfaces.set
                        (
                            intI,
                            LduGAMGInterface::New
                            (
                                fineMeshInterfaces[intI].type(),
                                intI,
                                dummyMeshInterfaces,
                                is
                            )
                        );
                    }
                }

                forAll(dummyPrimMeshInterfaces, intI)
                {
                    if (dummyPrimMeshInterfaces.set(intI))
                    {
                        dummyMeshInterfaces.set
                        (
                            intI,
                            &dummyPrimMeshInterfaces[intI]
                        );
                    }
                }

                // - pass in the incorrect (fine) mesh instead of the coarse
                // - pass in the dummy interfaces
                agglomerateMatrix
                (
                    fineLevelIndex,
                    agglomeration_.meshLevel(fineLevelIndex),
                    dummyMeshInterfaces
                );

                procAgglomerateMatrix
                (
                    agglomeration_.procAgglomMap(fineLevelIndex + 1),
                    agglomeration_.agglomProcIDs(fineLevelIndex + 1),
                    fineLevelIndex
                );
            }
            else
            {
                agglomerateMatrix
                (
                    fineLevelIndex,
                    agglomeration_.meshLevel(fineLevelIndex + 1),
                    agglomeration_.interfaceLevel(fineLevelIndex + 1)
                );
            }
        }
    }
    else
    {
        forAll(agglomeration_, fineLevelIndex)
        {
            agglomerateMatrix
            (
                fineLevelIndex,
                agglomeration_.meshLevel(fineLevelIndex + 1),
                agglomeration_.interfaceLevel(fineLevelIndex + 1)
            );
        }
    }

    if
    (
        (this->log_ >= 2)
     || (LduMatrix<Type, DType, LUType>::debug >= 2)
    )
    {
        for (label leveli = 0; leveli <= matrixLevels_.size(); ++leveli)
        {
            if (leveli == 0 || matrixLevels_.set(leveli - 1))
            {
                const auto& matrixLevelRef = matrixLevel(leveli);
                const auto& interfaces = interfaceLevel(leveli);

                Pout<< "level:" << leveli << nl
                    << "    nCells:" << matrixLevelRef.diag().size() << nl
                    << "    nFaces:" << matrixLevelRef.upper().size() << nl
                    << "    nInterfaces:" << interfaces.size()
                    << endl;
            }
            else
            {
                Pout<< "level:" << leveli << " : no matrix" << endl;
            }
        }

        Pout<< endl;
    }

    if (!matrixLevels_.size())
    {
        FatalErrorInFunction
            << "No coarse levels created, either matrix too small for GAMG"
               " or nCellsInCoarsestLevel too large.\n"
               "    Either choose another solver or reduce "
               "nCellsInCoarsestLevel."
            << exit(FatalError);
    }

    const label coarsestLevel = matrixLevels_.size() - 1;

    // On processor-agglomeration slaves the coarsest level is dormant
    // (null matrix); only the master constructs the coarsest-level solver.
    if (matrixLevels_.set(coarsestLevel))
    {
        if (directSolveCoarsest_)
        {
            buildCoarsestDirectMatrix();
            return;
        }

        dictionary coarseSolverControls;

        entry* coarseEntry = this->controlDict_.findEntry
        (
            "coarsestLevelCorr",
            keyType::LITERAL_RECURSIVE
        );

        if (coarseEntry && coarseEntry->isDict())
        {
            coarseSolverControls = coarseEntry->dict();
        }
        else
        {
            coarseSolverControls.set("solver", word("SmoothSolver"));
            label nCoarsestSweeps(nFinestSweeps_ > 0 ? nFinestSweeps_ : 1);

            if constexpr(!std::is_same<DType, scalar>::value)
            {
                nCoarsestSweeps = -nCoarsestSweeps;
            }

            coarseSolverControls.set
            (
                "nSweeps",
                nCoarsestSweeps
            );
            coarseSolverControls.set("tolerance", this->tolerance_);
            coarseSolverControls.set("relTol", this->relTol_);
            coarseSolverControls.set("minIter", label(1));
            coarseSolverControls.set("maxIter", this->maxIter_);
            coarseSolverControls.set
            (
                "smoother",
                this->controlDict_.template getOrDefault<word>
                (
                    "smoother",
                    "GaussSeidel"
                )
            );
        }

        if (coarseSolverControls.get<word>("solver") == typeName)
        {
            FatalIOErrorInFunction(coarseSolverControls)
                << "coarsestLevelCorr cannot select GAMG recursively"
                << exit(FatalIOError);
        }

        coarsestSolverPtr_ = LduMatrix<Type, DType, LUType>::solver::New
        (
            "coarsestLevelCorr",
            matrixLevels_[coarsestLevel],
            coarseSolverControls
        );
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::readControls()
{
    LduMatrix<Type, DType, LUType>::solver::readControls();

    this->controlDict_.readIfPresent("nPreSweeps", nPreSweeps_);
    this->controlDict_.readIfPresent
    (
        "preSweepsLevelMultiplier",
        preSweepsLevelMultiplier_
    );
    this->controlDict_.readIfPresent("maxPreSweeps", maxPreSweeps_);
    this->controlDict_.readIfPresent("nPostSweeps", nPostSweeps_);
    this->controlDict_.readIfPresent
    (
        "postSweepsLevelMultiplier",
        postSweepsLevelMultiplier_
    );
    this->controlDict_.readIfPresent("maxPostSweeps", maxPostSweeps_);
    this->controlDict_.readIfPresent("nFinestSweeps", nFinestSweeps_);
    this->controlDict_.readIfPresent
    (
        "interpolateCorrection",
        interpolateCorrection_
    );
    this->controlDict_.readIfPresent("scaleCorrection", scaleCorrection_);
    this->controlDict_.readIfPresent
    (
        "directSolveCoarsest",
        directSolveCoarsest_
    );

    if (this->log_ >= 2)
    {
        Info<< "LduGAMGSolver settings :"
            << " nPreSweeps:" << nPreSweeps_
            << " preSweepsLevelMultiplier:" << preSweepsLevelMultiplier_
            << " maxPreSweeps:" << maxPreSweeps_
            << " nPostSweeps:" << nPostSweeps_
            << " postSweepsLevelMultiplier:" << postSweepsLevelMultiplier_
            << " maxPostSweeps:" << maxPostSweeps_
            << " nFinestSweeps:" << nFinestSweeps_
            << " interpolateCorrection:" << interpolateCorrection_
            << " scaleCorrection:" << scaleCorrection_
            << " directSolveCoarsest:" << directSolveCoarsest_
            << endl;
    }
}


template<class Type, class DType, class LUType>
const Foam::LduMatrix<Type, DType, LUType>&
Foam::LduGAMGSolver<Type, DType, LUType>::matrixLevel
(
    const label leveli
) const
{
    return leveli ? matrixLevels_[leveli - 1] : this->matrix_;
}


template<class Type, class DType, class LUType>
const Foam::LduInterfaceFieldPtrsList<Type>&
Foam::LduGAMGSolver<Type, DType, LUType>::interfaceLevel
(
    const label leveli
) const
{
    return leveli ? interfaceLevels_[leveli - 1] : this->matrix_.interfaces();
}


template<class Type, class DType, class LUType>
const Foam::FieldField<Foam::Field, LUType>&
Foam::LduGAMGSolver<Type, DType, LUType>::interfaceUpperCoeffsLevel
(
    const label leveli
) const
{
    return
        leveli
      ? interfaceLevelsUpper_[leveli - 1]
      : this->matrix_.interfacesUpper();
}


template<class Type, class DType, class LUType>
const Foam::FieldField<Foam::Field, LUType>&
Foam::LduGAMGSolver<Type, DType, LUType>::interfaceLowerCoeffsLevel
(
    const label leveli
) const
{
    return
        leveli
      ? interfaceLevelsLower_[leveli - 1]
      : this->matrix_.interfacesLower();
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::agglomerateMatrix
(
    const label fineLevelIndex,
    const lduMesh& coarseMesh,
    const lduInterfacePtrsList& coarseMeshInterfaces
)
{
    const auto& fineMatrix = matrixLevel(fineLevelIndex);

    const label nCoarseFaces = agglomeration_.nFaces(fineLevelIndex);
    const label nCoarseCells = agglomeration_.nCells(fineLevelIndex);

    matrixLevels_.set
    (
        fineLevelIndex,
        new LduMatrix<Type, DType, LUType>(coarseMesh)
    );

    auto& coarseMatrix = matrixLevels_[fineLevelIndex];

    // Size with the cached coarse nCells and not the actual coarseMesh size
    // since the mesh may be a dummy when processor agglomerating.
    Field<DType>& coarseDiag = coarseMatrix.diag();
    coarseDiag.resize_nocopy(nCoarseCells);
    coarseDiag = Foam::zero{};

    agglomeration_.restrictField(coarseDiag, fineMatrix.diag(), fineLevelIndex);

    if (fineMatrix.hasSource())
    {
        agglomeration_.restrictField
        (
            coarseMatrix.source(),
            fineMatrix.source(),
            fineLevelIndex
        );
    }

    const auto& fineInterfaces = interfaceLevel(fineLevelIndex);

    primitiveInterfaceLevels_.set
    (
        fineLevelIndex,
        new PtrList<LduInterfaceField<Type>>(fineInterfaces.size())
    );

    interfaceLevels_.set
    (
        fineLevelIndex,
        new LduInterfaceFieldPtrsList<Type>(fineInterfaces.size())
    );

    interfaceLevelsUpper_.set
    (
        fineLevelIndex,
        new FieldField<Field, LUType>(fineInterfaces.size())
    );

    interfaceLevelsLower_.set
    (
        fineLevelIndex,
        new FieldField<Field, LUType>(fineInterfaces.size())
    );

    agglomerateInterfaceCoefficients
    (
        fineLevelIndex,
        coarseMeshInterfaces,
        primitiveInterfaceLevels_[fineLevelIndex],
        interfaceLevels_[fineLevelIndex],
        interfaceLevelsUpper_[fineLevelIndex],
        interfaceLevelsLower_[fineLevelIndex]
    );

    coarseMatrix.interfaces() = interfaceLevels_[fineLevelIndex];
    coarseMatrix.interfacesUpper() =
        FieldField<Field, LUType>(interfaceLevelsUpper_[fineLevelIndex]);
    coarseMatrix.interfacesLower() =
        FieldField<Field, LUType>(interfaceLevelsLower_[fineLevelIndex]);

    const labelList& faceRestrictAddr =
        agglomeration_.faceRestrictAddressing(fineLevelIndex);
    const boolList& faceFlipMap = agglomeration_.faceFlipMap(fineLevelIndex);

    if constexpr (std::is_same<DType, LUType>::value)
    {
        if (faceRestrictAddr.usePool() && faceFlipMap.usePool())
        {
            const bool hasLower = fineMatrix.hasLower();

            Field<LUType>& coarseUpper = coarseMatrix.upper();
            coarseUpper.resize_nocopy(nCoarseFaces);
            coarseUpper = Foam::zero{};

            foamExecutor exec;
            DType* __restrict__ coarseDiagPtr = coarseDiag.begin();
            LUType* __restrict__ coarseUpperPtr = coarseUpper.begin();
            const LUType* const __restrict__ fineUpperPtr =
                fineMatrix.upper().cbegin();
            const label* const __restrict__ faceRestrictPtr =
                faceRestrictAddr.cbegin();
            const bool* const __restrict__ faceFlipPtr = faceFlipMap.cbegin();

            if (hasLower)
            {
                Field<LUType>& coarseLower = coarseMatrix.lower();
                coarseLower.resize_nocopy(nCoarseFaces);
                coarseLower = Foam::zero{};

                LUType* __restrict__ coarseLowerPtr = coarseLower.begin();
                const LUType* const __restrict__ fineLowerPtr =
                    fineMatrix.lower().cbegin();

                auto Lambda = [=](label fineFacei)
                {
                    const label coarseFace = faceRestrictPtr[fineFacei];

                    if (coarseFace >= 0)
                    {
                        const bool flip = faceFlipPtr[fineFacei];

                        foamAtomic::AtomicAdd
                        (
                            coarseUpperPtr[coarseFace],
                            flip ? fineLowerPtr[fineFacei]
                                 : fineUpperPtr[fineFacei]
                        );
                        foamAtomic::AtomicAdd
                        (
                            coarseLowerPtr[coarseFace],
                            flip ? fineUpperPtr[fineFacei]
                                 : fineLowerPtr[fineFacei]
                        );
                    }
                    else
                    {
                        foamAtomic::AtomicAdd
                        (
                            coarseDiagPtr[-1 - coarseFace],
                            fineUpperPtr[fineFacei] + fineLowerPtr[fineFacei]
                        );
                    }
                };
                exec.parallelFor(Lambda, faceRestrictAddr.size());
            }
            else
            {
                auto Lambda = [=](label fineFacei)
                {
                    const label coarseFace = faceRestrictPtr[fineFacei];

                    if (coarseFace >= 0)
                    {
                        foamAtomic::AtomicAdd
                        (
                            coarseUpperPtr[coarseFace],
                            fineUpperPtr[fineFacei]
                        );
                    }
                    else
                    {
                        foamAtomic::AtomicAdd
                        (
                            coarseDiagPtr[-1 - coarseFace],
                            2*fineUpperPtr[fineFacei]
                        );
                    }
                };
                exec.parallelFor(Lambda, faceRestrictAddr.size());
            }

            return;
        }
    }

    if (fineMatrix.hasLower())
    {
        const Field<LUType>& fineUpper = fineMatrix.upper();
        const Field<LUType>& fineLower = fineMatrix.lower();

        Field<LUType>& coarseUpper = coarseMatrix.upper();
        Field<LUType>& coarseLower = coarseMatrix.lower();

        coarseUpper.resize_nocopy(nCoarseFaces);
        coarseLower.resize_nocopy(nCoarseFaces);
        coarseUpper = Foam::zero{};
        coarseLower = Foam::zero{};

        forAll(faceRestrictAddr, fineFacei)
        {
            const label coarseFace = faceRestrictAddr[fineFacei];

            if (coarseFace >= 0)
            {
                if (!faceFlipMap[fineFacei])
                {
                    coarseUpper[coarseFace] += fineUpper[fineFacei];
                    coarseLower[coarseFace] += fineLower[fineFacei];
                }
                else
                {
                    coarseUpper[coarseFace] += fineLower[fineFacei];
                    coarseLower[coarseFace] += fineUpper[fineFacei];
                }
            }
            else
            {
                if constexpr(std::is_same<DType, LUType>::value)
                {
                    coarseDiag[-1 - coarseFace] +=
                        fineUpper[fineFacei] + fineLower[fineFacei];
                }
                else
                {
                    DType diagContribution(pTraits<DType>::one);
                    diagContribution *=
                        fineUpper[fineFacei] + fineLower[fineFacei];
                    coarseDiag[-1 - coarseFace] += diagContribution;
                }
            }
        }
    }
    else
    {
        const Field<LUType>& fineUpper = fineMatrix.upper();
        Field<LUType>& coarseUpper = coarseMatrix.upper();

        coarseUpper.resize_nocopy(nCoarseFaces);
        coarseUpper = Foam::zero{};

        forAll(faceRestrictAddr, fineFacei)
        {
            const label coarseFace = faceRestrictAddr[fineFacei];

            if (coarseFace >= 0)
            {
                coarseUpper[coarseFace] += fineUpper[fineFacei];
            }
            else
            {
                if constexpr(std::is_same<DType, LUType>::value)
                {
                    coarseDiag[-1 - coarseFace] += 2*fineUpper[fineFacei];
                }
                else
                {
                    DType diagContribution(pTraits<DType>::one);
                    diagContribution *= 2*fineUpper[fineFacei];
                    coarseDiag[-1 - coarseFace] += diagContribution;
                }
            }
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::agglomerateInterfaceCoefficients
(
    const label fineLevelIndex,
    const lduInterfacePtrsList& coarseMeshInterfaces,
    PtrList<LduInterfaceField<Type>>& coarsePrimInterfaces,
    LduInterfaceFieldPtrsList<Type>& coarseInterfaces,
    FieldField<Field, LUType>& coarseInterfaceUpper,
    FieldField<Field, LUType>& coarseInterfaceLower
) const
{
    const auto& fineInterfaces = interfaceLevel(fineLevelIndex);
    const auto& fineInterfaceUpper = interfaceUpperCoeffsLevel(fineLevelIndex);
    const auto& fineInterfaceLower = interfaceLowerCoeffsLevel(fineLevelIndex);

    const labelListList& patchFineToCoarse =
        agglomeration_.patchFaceRestrictAddressing(fineLevelIndex);
    const labelList& nPatchFaces = agglomeration_.nPatchFaces(fineLevelIndex);

    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            const auto& coarseInterface =
                refCast<const LduGAMGInterface>(coarseMeshInterfaces[inti]);

            coarsePrimInterfaces.set
            (
                inti,
                newLduGAMGInterfaceField<Type>
                (
                    coarseInterface,
                    fineInterfaces[inti]
                ).ptr()
            );

            coarseInterfaces.set(inti, &coarsePrimInterfaces[inti]);

            const labelList& faceRestrictAddressing = patchFineToCoarse[inti];

            coarseInterfaceUpper.set
            (
                inti,
                new Field<LUType>(nPatchFaces[inti], Zero)
            );
            agglomeration_.restrictField
            (
                coarseInterfaceUpper[inti],
                fineInterfaceUpper[inti],
                faceRestrictAddressing
            );

            coarseInterfaceLower.set
            (
                inti,
                new Field<LUType>(nPatchFaces[inti], Zero)
            );

            if (fineInterfaceLower.set(inti))
            {
                agglomeration_.restrictField
                (
                    coarseInterfaceLower[inti],
                    fineInterfaceLower[inti],
                    faceRestrictAddressing
                );
            }
            else
            {
                agglomeration_.restrictField
                (
                    coarseInterfaceLower[inti],
                    fineInterfaceUpper[inti],
                    faceRestrictAddressing
                );
            }
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::gatherMatrices
(
    const label destLevel,
    const label comm,

    const LduMatrix<Type, DType, LUType>& mat,
    const FieldField<Field, LUType>& interfaceBouCoeffs,
    const FieldField<Field, LUType>& interfaceIntCoeffs,
    const LduInterfaceFieldPtrsList<Type>& interfaces,

    PtrList<LduMatrix<Type, DType, LUType>>& otherMats,
    PtrList<FieldField<Field, LUType>>& otherBouCoeffs,
    PtrList<FieldField<Field, LUType>>& otherIntCoeffs
) const
{
    PstreamBuffers pBufs(comm);

    // Send to master
    if (!UPstream::master(comm))
    {
        // All interfaceBouCoeffs need to be sent across
        bitSet validCoeffs(interfaces.size());
        forAll(interfaceBouCoeffs, intI)
        {
            if (interfaceBouCoeffs.set(intI))
            {
                validCoeffs.set(intI);
            }
        }

        UOPstream toMaster(UPstream::masterNo(), pBufs);

        // The templated LduMatrix has no Istream constructor, so the
        // individual fields are serialised explicitly.  The source is not
        // transferred: coarse-level sources are (re)set during the V-cycle.
        const bool hasDiag = mat.hasDiag();
        const bool hasUpper = mat.hasUpper();
        const bool hasLower = mat.hasLower();

        toMaster
            << hasDiag << token::SPACE
            << hasUpper << token::SPACE
            << hasLower;

        if (hasDiag)
        {
            toMaster << mat.diag();
        }
        if (hasUpper)
        {
            toMaster << mat.upper();
        }
        if (hasLower)
        {
            toMaster << mat.lower();
        }

        toMaster << validCoeffs;

        for (const label intI : validCoeffs)
        {
            toMaster
                << interfaceBouCoeffs[intI]
                << interfaceIntCoeffs[intI];
        }
    }

    // Wait for finish
    pBufs.finishedGathers();

    // Consume
    if (UPstream::master(comm))
    {
        const label nProcs = UPstream::nProcs(comm);

        const lduMesh& destMesh = agglomeration_.meshLevel(destLevel);

        otherMats.resize(nProcs - 1);
        otherBouCoeffs.resize(nProcs - 1);
        otherIntCoeffs.resize(nProcs - 1);

        for (const int proci : UPstream::subProcs(comm))
        {
            const label otherI = proci - 1;

            UIPstream fromProc(proci, pBufs);

            otherMats.set
            (
                otherI,
                new LduMatrix<Type, DType, LUType>(destMesh)
            );
            auto& otherMat = otherMats[otherI];

            bool hasDiag = false;
            bool hasUpper = false;
            bool hasLower = false;

            fromProc >> hasDiag >> hasUpper >> hasLower;

            if (hasDiag)
            {
                fromProc >> otherMat.diag();
            }
            if (hasUpper)
            {
                fromProc >> otherMat.upper();
            }
            if (hasLower)
            {
                fromProc >> otherMat.lower();
            }

            // Receive bit-set of valid interfaceCoeffs
            const bitSet validCoeffs(fromProc);

            otherBouCoeffs.emplace_set(otherI, validCoeffs.size());
            otherIntCoeffs.emplace_set(otherI, validCoeffs.size());

            // Receive individual interface contributions
            for (const label intI : validCoeffs)
            {
                otherBouCoeffs[otherI].emplace_set(intI, fromProc);
                otherIntCoeffs[otherI].emplace_set(intI, fromProc);
            }
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::procAgglomerateMatrix
(
    const labelList& procAgglomMap,
    const labelList& agglomProcIDs,
    const label levelI,

    autoPtr<LduMatrix<Type, DType, LUType>>& allMatrixPtr,
    FieldField<Field, LUType>& allInterfaceBouCoeffs,
    FieldField<Field, LUType>& allInterfaceIntCoeffs,
    PtrList<LduInterfaceField<Type>>& allPrimitiveInterfaces,
    LduInterfaceFieldPtrsList<Type>& allInterfaces
) const
{
    const LduMatrix<Type, DType, LUType>& coarsestMatrix =
        matrixLevels_[levelI];
    const LduInterfaceFieldPtrsList<Type>& coarsestInterfaces =
        interfaceLevels_[levelI];
    const FieldField<Field, LUType>& coarsestBouCoeffs =
        interfaceLevelsUpper_[levelI];
    const FieldField<Field, LUType>& coarsestIntCoeffs =
        interfaceLevelsLower_[levelI];

    // Communicator containing all processors to combine (=agglomProcIDs).
    // Result will be on master of communicator.
    const label agglomComm = agglomeration_.agglomCommunicator(levelI + 1);

    // Gather all matrix coefficients onto agglomProcIDs[0]
    PtrList<LduMatrix<Type, DType, LUType>> otherMats;
    PtrList<FieldField<Field, LUType>> otherBouCoeffs;
    PtrList<FieldField<Field, LUType>> otherIntCoeffs;
    gatherMatrices
    (
        levelI + 1,             // allMesh level (only on master)
        agglomComm,

        coarsestMatrix,
        coarsestBouCoeffs,
        coarsestIntCoeffs,
        coarsestInterfaces,

        otherMats,
        otherBouCoeffs,
        otherIntCoeffs
    );

    if (UPstream::master(agglomComm))
    {
        const lduMesh& allMesh = agglomeration_.meshLevel(levelI + 1);
        const labelList& cellOffsets = agglomeration_.cellOffsets(levelI + 1);
        const labelListList& faceMap = agglomeration_.faceMap(levelI + 1);
        const labelListList& boundaryMap =
            agglomeration_.boundaryMap(levelI + 1);
        const labelListListList& boundaryFaceMap =
            agglomeration_.boundaryFaceMap(levelI + 1);

        allMatrixPtr.reset(new LduMatrix<Type, DType, LUType>(allMesh));
        auto& allMatrix = allMatrixPtr();

        // Assemble matrix diagonal, source and off-diagonal coefficients
        if (coarsestMatrix.hasDiag())
        {
            Field<DType>& allDiag = allMatrix.diag();

            SubList<DType>(allDiag, coarsestMatrix.diag().size()) =
                coarsestMatrix.diag();

            forAll(otherMats, i)
            {
                SubList<DType>
                (
                    allDiag,
                    otherMats[i].diag().size(),
                    cellOffsets[i + 1]
                ) = otherMats[i].diag();
            }
        }
        if (coarsestMatrix.hasLower())
        {
            Field<LUType>& allLower = allMatrix.lower();
            UIndirectList<LUType>(allLower, faceMap[0]) =
                coarsestMatrix.lower();
            forAll(otherMats, i)
            {
                UIndirectList<LUType>(allLower, faceMap[i + 1]) =
                    otherMats[i].lower();
            }
        }
        if (coarsestMatrix.hasUpper())
        {
            Field<LUType>& allUpper = allMatrix.upper();
            UIndirectList<LUType>(allUpper, faceMap[0]) =
                coarsestMatrix.upper();
            forAll(otherMats, i)
            {
                UIndirectList<LUType>(allUpper, faceMap[i + 1]) =
                    otherMats[i].upper();
            }
        }

        // Agglomerate interface fields and coefficients
        lduInterfacePtrsList allMeshInterfaces = allMesh.interfaces();

        allInterfaceBouCoeffs.setSize(allMeshInterfaces.size());
        allInterfaceIntCoeffs.setSize(allMeshInterfaces.size());
        allPrimitiveInterfaces.setSize(allMeshInterfaces.size());
        allInterfaces.setSize(allMeshInterfaces.size());

        forAll(allMeshInterfaces, intI)
        {
            const lduInterface& patch = allMeshInterfaces[intI];
            const label size = patch.faceCells().size();

            allInterfaceBouCoeffs.set(intI, new Field<LUType>(size, Zero));
            allInterfaceIntCoeffs.set(intI, new Field<LUType>(size, Zero));
        }

        forAll(boundaryMap, proci)
        {
            const FieldField<Field, LUType>& procBouCoeffs
            (
                (proci == 0) ? coarsestBouCoeffs : otherBouCoeffs[proci - 1]
            );
            const FieldField<Field, LUType>& procIntCoeffs
            (
                (proci == 0) ? coarsestIntCoeffs : otherIntCoeffs[proci - 1]
            );

            const labelList& bMap = boundaryMap[proci];
            forAll(bMap, procIntI)
            {
                const label allIntI = bMap[procIntI];

                if (allIntI != -1)
                {
                    // This boundary has been preserved
                    if (!allInterfaces.set(allIntI))
                    {
                        if (!isA<processorLduGAMGInterface>
                            (allMeshInterfaces[allIntI]))
                        {
                            FatalErrorInFunction
                                << "Processor agglomeration of templated"
                                << " LduMatrix GAMG only supports surviving"
                                << " processor interfaces.  Interface "
                                << allIntI << " has type "
                                << allMeshInterfaces[allIntI].type()
                                << exit(FatalError);
                        }

                        const auto& intf =
                            refCast<const processorLduGAMGInterface>
                            (
                                allMeshInterfaces[allIntI]
                            );

                        // Recreate a fresh processor interface field.  All
                        // surviving interfaces are processor ones, hence no
                        // transform is required.
                        allPrimitiveInterfaces.set
                        (
                            allIntI,
                            new ProcessorLduGAMGInterfaceField<Type>
                            (
                                intf,
                                false,                  // doTransform
                                pTraits<Type>::rank     // rank
                            )
                        );

                        allInterfaces.set
                        (
                            allIntI,
                            &allPrimitiveInterfaces[allIntI]
                        );
                    }

                    // Map data from processor to complete mesh
                    Field<LUType>& allBou = allInterfaceBouCoeffs[allIntI];
                    Field<LUType>& allInt = allInterfaceIntCoeffs[allIntI];

                    const labelList& map = boundaryFaceMap[proci][procIntI];

                    const Field<LUType>& procBou = procBouCoeffs[procIntI];
                    const Field<LUType>& procInt = procIntCoeffs[procIntI];

                    forAll(map, i)
                    {
                        const label allFacei = map[i];
                        if (allFacei < 0)
                        {
                            FatalErrorInFunction
                                << "problem." << abort(FatalError);
                        }
                        allBou[allFacei] = procBou[i];
                        allInt[allFacei] = procInt[i];
                    }
                }
                else if (procBouCoeffs.set(procIntI))
                {
                    // Boundary has become an internal face
                    const labelList& map = boundaryFaceMap[proci][procIntI];
                    const Field<LUType>& procBou = procBouCoeffs[procIntI];
                    const Field<LUType>& procInt = procIntCoeffs[procIntI];

                    forAll(map, i)
                    {
                        if (map[i] >= 0)
                        {
                            const label allFacei = map[i];

                            if (coarsestMatrix.hasUpper())
                            {
                                allMatrix.upper()[allFacei] = -procBou[i];
                            }
                            if (coarsestMatrix.hasLower())
                            {
                                allMatrix.lower()[allFacei] = -procInt[i];
                            }
                        }
                        else
                        {
                            const label allFacei = -map[i] - 1;

                            if (coarsestMatrix.hasUpper())
                            {
                                allMatrix.upper()[allFacei] = -procInt[i];
                            }
                            if (coarsestMatrix.hasLower())
                            {
                                allMatrix.lower()[allFacei] = -procBou[i];
                            }
                        }
                    }
                }
            }
        }

        allMatrix.interfaces() = allInterfaces;
        allMatrix.interfacesUpper() =
            FieldField<Field, LUType>(allInterfaceBouCoeffs);
        allMatrix.interfacesLower() =
            FieldField<Field, LUType>(allInterfaceIntCoeffs);
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::procAgglomerateMatrix
(
    const labelList& procAgglomMap,
    const labelList& agglomProcIDs,
    const label levelI
)
{
    autoPtr<LduMatrix<Type, DType, LUType>> allMatrixPtr;
    autoPtr<FieldField<Field, LUType>> allInterfaceBouCoeffs
    (
        new FieldField<Field, LUType>(0)
    );
    autoPtr<FieldField<Field, LUType>> allInterfaceIntCoeffs
    (
        new FieldField<Field, LUType>(0)
    );
    autoPtr<PtrList<LduInterfaceField<Type>>> allPrimitiveInterfaces
    (
        new PtrList<LduInterfaceField<Type>>(0)
    );
    autoPtr<LduInterfaceFieldPtrsList<Type>> allInterfaces
    (
        new LduInterfaceFieldPtrsList<Type>(0)
    );

    procAgglomerateMatrix
    (
        procAgglomMap,
        agglomProcIDs,
        levelI,

        allMatrixPtr,
        allInterfaceBouCoeffs(),
        allInterfaceIntCoeffs(),
        allPrimitiveInterfaces(),
        allInterfaces()
    );

    matrixLevels_.set(levelI, allMatrixPtr);
    interfaceLevelsUpper_.set(levelI, allInterfaceBouCoeffs);
    interfaceLevelsLower_.set(levelI, allInterfaceIntCoeffs);
    primitiveInterfaceLevels_.set(levelI, allPrimitiveInterfaces);
    interfaceLevels_.set(levelI, allInterfaces);
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::scale
(
    Field<Type>& field,
    Field<Type>& Acf,
    const LduMatrix<Type, DType, LUType>& A,
    const Field<Type>& source
) const
{
    A.Amul(Acf, field);

    const label comm = A.mesh().comm();

    const scalar sf =
    (
        gSumProd(field, source, comm)
      / stabilise(gSumProd(field, Acf, comm), pTraits<scalar>::vsmall)
    );

    if (LduMatrix<Type, DType, LUType>::debug >= 2)
    {
        Info<< sf << ' ';
    }

    Type* __restrict__ fieldPtr = field.begin();
    const Type* const __restrict__ sourcePtr = source.begin();
    const Type* const __restrict__ AcfPtr = Acf.begin();
    const DType* const __restrict__ diagPtr = A.diag().begin();

    foamExecutor exec;
    const label nCells = field.size();

    auto LambdaScale = [=](label celli)
    {
        const Type residual = sourcePtr[celli] - sf*AcfPtr[celli];

        if constexpr(std::is_same<DType, Type>::value)
        {
            fieldPtr[celli] =
                sf*fieldPtr[celli] + cmptDivide(residual, diagPtr[celli]);
        }
        else
        {
            if constexpr
            (
                std::is_same<Type, vector4>::value
             && std::is_same<DType, tensor4>::value
            )
            {
                fieldPtr[celli] =
                    sf*fieldPtr[celli]
                  + ::Foam::solve(diagPtr[celli], residual);
            }
            else
            {
                fieldPtr[celli] = sf*fieldPtr[celli] + residual/diagPtr[celli];
            }
        }
    };
    exec.parallelFor(LambdaScale, nCells);
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::interpolate
(
    Field<Type>& psi,
    Field<Type>& Apsi,
    const LduMatrix<Type, DType, LUType>& A
) const
{
    Type* __restrict__ psiPtr = psi.begin();
    Type* __restrict__ ApsiPtr = Apsi.begin();

    const DType* const __restrict__ diagPtr = A.diag().begin();
    const LUType* const __restrict__ upperPtr = A.upper().begin();
    const LUType* const __restrict__ lowerPtr = A.lower().begin();
    const label* const __restrict__ uPtr = A.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = A.lduAddr().lowerAddr().begin();

    Apsi = Zero;

    const label startRequest = UPstream::nRequests();

    A.initMatrixInterfaces
    (
        true,
        A.interfacesUpper(),
        psi,
        Apsi
    );

    foamExecutor exec;
    const label nFaces = A.upper().size();

    auto LambdaOffDiag = [=](label facei)
    {
        foamAtomic::AtomicAdd
        (
            ApsiPtr[uPtr[facei]],
            dot(lowerPtr[facei], psiPtr[lPtr[facei]])
        );
        foamAtomic::AtomicAdd
        (
            ApsiPtr[lPtr[facei]],
            dot(upperPtr[facei], psiPtr[uPtr[facei]])
        );
    };
    exec.parallelFor(LambdaOffDiag, nFaces);

    A.updateMatrixInterfaces
    (
        true,
        A.interfacesUpper(),
        psi,
        Apsi,
        startRequest
    );

    const label nCells = psi.size();

    auto LambdaDiag = [=](label celli)
    {
        if constexpr(std::is_same<DType, Type>::value)
        {
            psiPtr[celli] = -cmptDivide(ApsiPtr[celli], diagPtr[celli]);
        }
        else
        {
            if constexpr
            (
                std::is_same<Type, vector4>::value
             && std::is_same<DType, tensor4>::value
            )
            {
                psiPtr[celli] =
                    -::Foam::solve(diagPtr[celli], ApsiPtr[celli]);
            }
            else
            {
                psiPtr[celli] = -ApsiPtr[celli]/diagPtr[celli];
            }
        }
    };
    exec.parallelFor(LambdaDiag, nCells);
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::interpolate
(
    Field<Type>& psi,
    Field<Type>& Apsi,
    const LduMatrix<Type, DType, LUType>& A,
    const labelList& restrictAddressing,
    const Field<Type>& psiC
) const
{
    interpolate(psi, Apsi, A);

    Type* __restrict__ psiPtr = psi.begin();
    const DType* const __restrict__ diagPtr = A.diag().begin();
    const Type* const __restrict__ psiCPtr = psiC.begin();

    Field<Type> corrC(psiC.size(), Zero);
    Field<DType> diagC(psiC.size(), Zero);

    Type* __restrict__ corrCPtr = corrC.begin();
    DType* __restrict__ diagCPtr = diagC.begin();
    const label* const __restrict__ restrictAddressingPtr =
        restrictAddressing.begin();

    foamExecutor exec;
    const label nCells = psi.size();

    auto LambdaRestrict = [=](label celli)
    {
        const label coarseCelli = restrictAddressingPtr[celli];

        Type diagPsi;
        if constexpr(std::is_same<DType, Type>::value)
        {
            diagPsi = cmptMultiply(diagPtr[celli], psiPtr[celli]);
        }
        else
        {
            diagPsi = dot(diagPtr[celli], psiPtr[celli]);
        }

        foamAtomic::AtomicAdd
        (
            corrCPtr[coarseCelli],
            diagPsi
        );
        foamAtomic::AtomicAdd(diagCPtr[coarseCelli], diagPtr[celli]);
    };
    exec.parallelFor(LambdaRestrict, nCells);

    const label nCoarseCells = corrC.size();

    auto LambdaCoarse = [=](label coarseCelli)
    {
        if constexpr(std::is_same<DType, Type>::value)
        {
            corrCPtr[coarseCelli] =
                psiCPtr[coarseCelli]
              - cmptDivide(corrCPtr[coarseCelli], diagCPtr[coarseCelli]);
        }
        else
        {
            if constexpr
            (
                std::is_same<Type, vector4>::value
             && std::is_same<DType, tensor4>::value
            )
            {
                corrCPtr[coarseCelli] =
                    psiCPtr[coarseCelli]
                  - ::Foam::solve(diagCPtr[coarseCelli], corrCPtr[coarseCelli]);
            }
            else
            {
                corrCPtr[coarseCelli] =
                    psiCPtr[coarseCelli]
                  - corrCPtr[coarseCelli]/diagCPtr[coarseCelli];
            }
        }
    };
    exec.parallelFor(LambdaCoarse, nCoarseCells);

    auto LambdaProlong = [=](label celli)
    {
        psiPtr[celli] += corrCPtr[restrictAddressingPtr[celli]];
    };
    exec.parallelFor(LambdaProlong, nCells);
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::initVcycle
(
    PtrList<Field<Type>>& coarseCorrFields,
    PtrList<Field<Type>>& coarseSources,
    PtrList<typename LduMatrix<Type, DType, LUType>::smoother>& smoothers
) const
{
    dictionary smootherControls(this->controlDict_);
    smootherControls.set
    (
        "smoother",
        this->controlDict_.template getOrDefault<word>
        (
            "smoother",
            "GaussSeidel"
        )
    );

    coarseCorrFields.setSize(matrixLevels_.size());
    coarseSources.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    smoothers.set
    (
        0,
        LduMatrix<Type, DType, LUType>::smoother::New
        (
            this->fieldName_,
            this->matrix_,
            smootherControls
        )
    );

    forAll(matrixLevels_, leveli)
    {
        // Source fields exist on every processor that contributes to a level
        // (including processor-agglomeration slaves whose matrix is dormant)
        // so that the gather/scatter collectives are matched.
        if (agglomeration_.nCells(leveli) >= 0)
        {
            coarseSources.set
            (
                leveli,
                new Field<Type>(agglomeration_.nCells(leveli), Zero)
            );
        }

        if (matrixLevels_.set(leveli))
        {
            const label nCoarseCells = matrixLevels_[leveli].diag().size();

            coarseCorrFields.set(leveli, new Field<Type>(nCoarseCells, Zero));

            smoothers.set
            (
                leveli + 1,
                LduMatrix<Type, DType, LUType>::smoother::New
                (
                    this->fieldName_,
                    matrixLevels_[leveli],
                    smootherControls
                )
            );
        }
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::Vcycle
(
    const PtrList<typename LduMatrix<Type, DType, LUType>::smoother>& smoothers,
    Field<Type>& psi,
    Field<Type>& Apsi,
    Field<Type>& finestCorrection,
    Field<Type>& finestResidual,
    PtrList<Field<Type>>& coarseCorrFields,
    PtrList<Field<Type>>& coarseSources
) const
{
    const label coarsestLevel = matrixLevels_.size() - 1;
    const label finestSweeps = (nFinestSweeps_ > 0 ? nFinestSweeps_ : 1);

    agglomeration_.restrictField(coarseSources[0], finestResidual, 0, true);

    if
    (
        nPreSweeps_
     && ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
    )
    {
        Info<< "Pre-smoothing scaling factors: ";
    }

    for (label leveli = 0; leveli < coarsestLevel; ++leveli)
    {
        if (coarseSources.set(leveli + 1))
        {
            auto& coarseMatrix = matrixLevels_[leveli];
            coarseMatrix.source() = coarseSources[leveli];

            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] = Zero;

                smoothers[leveli + 1].smooth
                (
                    coarseCorrFields[leveli],
                    min
                    (
                        nPreSweeps_ + preSweepsLevelMultiplier_*leveli,
                        maxPreSweeps_
                    )
                );

                if (scaleCorrection_ && leveli < coarsestLevel - 1)
                {
                    Field<Type> Acf(coarseCorrFields[leveli].size());

                    scale
                    (
                        coarseCorrFields[leveli],
                        Acf,
                        coarseMatrix,
                        coarseSources[leveli]
                    );
                }

                coarseMatrix.source() = coarseSources[leveli];
                coarseMatrix.residual
                (
                    coarseSources[leveli],
                    coarseCorrFields[leveli]
                );
            }

            agglomeration_.restrictField
            (
                coarseSources[leveli + 1],
                coarseSources[leveli],
                leveli + 1,
                true
            );
        }
    }

    if
    (
        nPreSweeps_
     && ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
    )
    {
        Info<< endl;
    }

    // Solve coarsest level (skipped on processor-agglomeration slaves where
    // the coarsest matrix is dormant)
    if (coarseCorrFields.set(coarsestLevel))
    {
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    if ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
    {
        Info<< "Post-smoothing scaling factors: ";
    }
    Field<Type> dummyField(0);
    Field<Type> work;

    for (label leveli = coarsestLevel - 1; leveli >= 0; --leveli)
    {
        if (coarseCorrFields.set(leveli))
        {
            Field<Type> preSmoothedCoarseCorrField;

            if (nPreSweeps_)
            {
                preSmoothedCoarseCorrField = coarseCorrFields[leveli];
            }

            const auto& cf = agglomeration_.prolongField
            (
                coarseCorrFields[leveli],
                work,
                (
                    coarseCorrFields.set(leveli + 1)
                  ? coarseCorrFields[leveli + 1]
                  : dummyField
                ),
                leveli + 1
            );

            Field<Type> Acf(coarseCorrFields[leveli].size(), Zero);

            if (interpolateCorrection_)
            {
                interpolate
                (
                    coarseCorrFields[leveli],
                    Acf,
                    matrixLevels_[leveli],
                    agglomeration_.restrictAddressing(leveli + 1),
                    cf
                );
            }

            if
            (
                scaleCorrection_
             && (interpolateCorrection_ || leveli < coarsestLevel - 1)
            )
            {
                scale
                (
                    coarseCorrFields[leveli],
                    Acf,
                    matrixLevels_[leveli],
                    coarseSources[leveli]
                );
            }

            if (nPreSweeps_)
            {
                coarseCorrFields[leveli] += preSmoothedCoarseCorrField;
            }

            matrixLevels_[leveli].source() = coarseSources[leveli];

            smoothers[leveli + 1].smooth
            (
                coarseCorrFields[leveli],
                min
                (
                    nPostSweeps_ + postSweepsLevelMultiplier_*leveli,
                    maxPostSweeps_
                )
            );
        }
    }

    if ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
    {
        Info<< endl;
    }

    agglomeration_.prolongField(finestCorrection, coarseCorrFields[0], 0, true);

    if (interpolateCorrection_)
    {
        interpolate
        (
            finestCorrection,
            Apsi,
            this->matrix_,
            agglomeration_.restrictAddressing(0),
            coarseCorrFields[0]
        );
    }

    if (scaleCorrection_)
    {
        scale(finestCorrection, Apsi, this->matrix_, finestResidual);
    }

    psi += finestCorrection;

    smoothers[0].smooth(psi, finestSweeps);
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::solveCoarsestLevel
(
    Field<Type>& coarsestCorrField,
    const Field<Type>& coarsestSource
) const
{
    if (directSolveCoarsest_)
    {
        directSolveCoarsestLevel(coarsestCorrField, coarsestSource);
        return;
    }

    const label coarsestLevel = matrixLevels_.size() - 1;

    matrixLevels_[coarsestLevel].source() = coarsestSource;
    coarsestCorrField = Zero;

    const SolverPerformance<Type> coarseSolverPerf
    (
        coarsestSolverPtr_->solve(coarsestCorrField)
    );

    if
    (
        (this->log_ >= 2)
     || LduMatrix<Type, DType, LUType>::debug
    )
    {
        coarseSolverPerf.print(Info);
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::buildCoarsestDirectMatrix()
{
    if constexpr
    (
        !std::is_same<DType, scalar>::value
     || !std::is_same<LUType, scalar>::value
    )
    {
        FatalErrorInFunction
            << "directSolveCoarsest is only supported for scalar "
            << "LduMatrix diagonal and off-diagonal coefficients"
            << exit(FatalError);

        return;
    }
    else
    {
    // Assemble the (possibly distributed) coarsest scalar matrix onto the
    // master of the coarsest mesh communicator and LU-factor it.  This is a
    // native equivalent of the legacy LUscalarMatrix assembly: the templated
    // GAMG matrix has scalar diagonal/off-diagonal coefficients (DType and
    // LUType are scalar), so the gathered system is a plain scalar dense
    // matrix.  The Type source is solved component-wise via LUBacksubstitute.
    const label coarsestLevel = matrixLevels_.size() - 1;

    const LduMatrix<Type, DType, LUType>& mat = matrixLevels_[coarsestLevel];
    const LduInterfaceFieldPtrsList<Type>& interfaces =
        interfaceLevels_[coarsestLevel];
    const FieldField<Field, LUType>& bouCoeffs =
        interfaceLevelsUpper_[coarsestLevel];

    // All coarse reductions/communications use the matrix mesh communicator
    coarsestComm_ = mat.mesh().comm();
    const label comm = coarsestComm_;

    // Local scalar matrix snapshot
    const labelList localUpperAddr(mat.lduAddr().upperAddr());
    const labelList localLowerAddr(mat.lduAddr().lowerAddr());
    const scalarField localDiag(mat.diag());
    scalarField localUpper;
    scalarField localLower;
    if (mat.hasUpper() || mat.hasLower())
    {
        localUpper = mat.upper();
        localLower = mat.lower();
    }

    // Local interface snapshot (preserving the full interface indexing so
    // that cyclic neighbPatchID references remain valid on the master)
    const label nIf = interfaces.size();
    labelList localIfType(nIf, -1);             // -1 invalid, 0 proc, 1 cyclic
    List<labelList> localIfFaceCells(nIf);
    List<scalarField> localIfCoeffs(nIf);
    labelList localIfMyProc(nIf, -1);
    labelList localIfNeiProc(nIf, -1);
    labelList localIfTag(nIf, -1);
    labelList localIfNbrPatch(nIf, -1);

    forAll(interfaces, i)
    {
        if (!interfaces.set(i) || !bouCoeffs.set(i))
        {
            continue;
        }

        const lduInterface& li = interfaces[i].interface();
        localIfFaceCells[i] = li.faceCells();
        localIfCoeffs[i] = bouCoeffs[i];

        if (isA<processorLduInterface>(li))
        {
            const processorLduInterface& pli =
                refCast<const processorLduInterface>(li);

            localIfType[i] = 0;
            localIfMyProc[i] = pli.myProcNo();
            localIfNeiProc[i] = pli.neighbProcNo();
            localIfTag[i] = pli.tag();
        }
        else if (isA<cyclicLduInterface>(li))
        {
            const cyclicLduInterface& cli =
                refCast<const cyclicLduInterface>(li);

            localIfType[i] = 1;
            localIfNbrPatch[i] = cli.neighbPatchID();
        }
    }

    const label nProcs = UPstream::nProcs(comm);
    const label myProc = UPstream::myProcNo(comm);

    // Gather every rank's snapshot onto the master
    PstreamBuffers pBufs(comm);

    if (!UPstream::master(comm))
    {
        UOPstream toMaster(UPstream::masterNo(), pBufs);
        toMaster
            << localUpperAddr << localLowerAddr
            << localDiag << localUpper << localLower
            << localIfType << localIfFaceCells << localIfCoeffs
            << localIfMyProc << localIfNeiProc << localIfTag
            << localIfNbrPatch;
    }

    pBufs.finishedGathers();

    if (!UPstream::master(comm))
    {
        // Non-master ranks only participate in the coarsest solve scatter
        return;
    }

    labelListList allUpperAddr(nProcs);
    labelListList allLowerAddr(nProcs);
    List<scalarField> allDiag(nProcs);
    List<scalarField> allUpper(nProcs);
    List<scalarField> allLower(nProcs);
    labelListList allIfType(nProcs);
    List<List<labelList>> allIfFaceCells(nProcs);
    List<List<scalarField>> allIfCoeffs(nProcs);
    labelListList allIfMyProc(nProcs);
    labelListList allIfNeiProc(nProcs);
    labelListList allIfTag(nProcs);
    labelListList allIfNbrPatch(nProcs);

    // Master's own contribution
    allUpperAddr[myProc] = localUpperAddr;
    allLowerAddr[myProc] = localLowerAddr;
    allDiag[myProc] = localDiag;
    allUpper[myProc] = localUpper;
    allLower[myProc] = localLower;
    allIfType[myProc] = localIfType;
    allIfFaceCells[myProc] = localIfFaceCells;
    allIfCoeffs[myProc] = localIfCoeffs;
    allIfMyProc[myProc] = localIfMyProc;
    allIfNeiProc[myProc] = localIfNeiProc;
    allIfTag[myProc] = localIfTag;
    allIfNbrPatch[myProc] = localIfNbrPatch;

    for (const int proci : UPstream::subProcs(comm))
    {
        UIPstream fromProc(proci, pBufs);
        fromProc
            >> allUpperAddr[proci] >> allLowerAddr[proci]
            >> allDiag[proci] >> allUpper[proci] >> allLower[proci]
            >> allIfType[proci] >> allIfFaceCells[proci] >> allIfCoeffs[proci]
            >> allIfMyProc[proci] >> allIfNeiProc[proci] >> allIfTag[proci]
            >> allIfNbrPatch[proci];
    }

    // Cell offsets per rank
    coarsestProcOffsets_.resize_nocopy(nProcs + 1);
    coarsestProcOffsets_[0] = 0;
    for (label proci = 0; proci < nProcs; ++proci)
    {
        coarsestProcOffsets_[proci + 1] =
            coarsestProcOffsets_[proci] + allDiag[proci].size();
    }
    const label nTotal = coarsestProcOffsets_[nProcs];

    coarsestLUPtr_.reset(new scalarSquareMatrix(nTotal, Foam::zero{}));
    scalarSquareMatrix& m = coarsestLUPtr_();

    for (label proci = 0; proci < nProcs; ++proci)
    {
        const label offset = coarsestProcOffsets_[proci];

        const labelList& uAddr = allUpperAddr[proci];
        const labelList& lAddr = allLowerAddr[proci];
        const scalarField& diag = allDiag[proci];
        const scalarField& upper = allUpper[proci];
        const scalarField& lower = allLower[proci];

        const label nCells = diag.size();
        const label nFaces = upper.size();

        for (label cell = 0; cell < nCells; ++cell)
        {
            m(cell + offset, cell + offset) = diag[cell];
        }

        for (label face = 0; face < nFaces; ++face)
        {
            const label uCell = uAddr[face] + offset;
            const label lCell = lAddr[face] + offset;

            m(uCell, lCell) = lower[face];
            m(lCell, uCell) = upper[face];
        }

        // Interface (boundary) coupling
        const labelList& ifType = allIfType[proci];
        const List<labelList>& ifFaceCells = allIfFaceCells[proci];
        const List<scalarField>& ifCoeffs = allIfCoeffs[proci];

        forAll(ifType, inti)
        {
            if (ifType[inti] == 1)
            {
                // Cyclic interface: couple to its neighbour patch on the
                // same rank using the neighbour patch coefficients
                const label nbrInt = allIfNbrPatch[proci][inti];

                const labelList& lCells = ifFaceCells[inti];
                const labelList& uCells = ifFaceCells[nbrInt];
                const scalarField& nbrCoeffs = ifCoeffs[nbrInt];

                forAll(lCells, face)
                {
                    const label uCell = lCells[face] + offset;
                    const label lCell = uCells[face] + offset;

                    m(uCell, lCell) -= nbrCoeffs[face];
                }
            }
            else if (ifType[inti] == 0)
            {
                // Processor interface: couple to the matching interface on
                // the neighbour rank (only assembled once, from the lower
                // ranked processor)
                const label myP = allIfMyProc[proci][inti];
                const label neiP = allIfNeiProc[proci][inti];

                if (myP >= neiP)
                {
                    continue;
                }

                // Find the matching interface on the neighbour rank
                const labelList& neiType = allIfType[neiP];
                const labelList& neiMyProc = allIfMyProc[neiP];
                const labelList& neiNeiProc = allIfNeiProc[neiP];
                const labelList& neiTag = allIfTag[neiP];
                const label myTag = allIfTag[proci][inti];

                label neiInt = -1;
                forAll(neiType, j)
                {
                    if
                    (
                        neiType[j] == 0
                     && neiNeiProc[j] == myP
                     && neiMyProc[j] == neiP
                     && neiTag[j] == myTag
                    )
                    {
                        neiInt = j;
                        break;
                    }
                }

                if (neiInt == -1)
                {
                    FatalErrorInFunction
                        << "Could not find matching processor interface for"
                        << " coarsest direct solve (proc " << myP
                        << " -> " << neiP << ", tag " << myTag << ")"
                        << exit(FatalError);
                }

                const labelList& uCells = ifFaceCells[inti];
                const labelList& lCells = allIfFaceCells[neiP][neiInt];
                const scalarField& upperPtr = ifCoeffs[inti];
                const scalarField& lowerPtr = allIfCoeffs[neiP][neiInt];

                const label neiOffset = coarsestProcOffsets_[neiP];

                forAll(uCells, face)
                {
                    const label uCell = uCells[face] + offset;
                    const label lCell = lCells[face] + neiOffset;

                    m(uCell, lCell) -= lowerPtr[face];
                    m(lCell, uCell) -= upperPtr[face];
                }
            }
        }
    }

    LUDecompose(m, coarsestPivots_);
    }
}


template<class Type, class DType, class LUType>
void Foam::LduGAMGSolver<Type, DType, LUType>::directSolveCoarsestLevel
(
    Field<Type>& coarsestCorrField,
    const Field<Type>& coarsestSource
) const
{
    const label comm = coarsestComm_;
    const auto tag = UPstream::msgType();

    // Initialise the working solution with the source
    coarsestCorrField = coarsestSource;
    Field<Type>& x = coarsestCorrField;

    if (UPstream::parRun())
    {
        List<Type> allx;  // scratch space (on master)

        const label startOfRequests = UPstream::nRequests();

        // Gather the source onto the master (globalIndex::gather pattern)
        if (UPstream::master(comm))
        {
            allx.resize(coarsestProcOffsets_.last());

            SubList<Type>(allx, x.size()) = x;

            for (const int proci : UPstream::subProcs(comm))
            {
                SubList<Type> procSlot
                (
                    allx,
                    coarsestProcOffsets_[proci+1]-coarsestProcOffsets_[proci],
                    coarsestProcOffsets_[proci]
                );

                if (procSlot.empty())
                {
                    // Nothing to do
                }
                else if (is_contiguous<Type>::value)
                {
                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        procSlot.data_bytes(),
                        procSlot.size_bytes(),
                        tag,
                        comm
                    );
                }
                else
                {
                    IPstream::recv(procSlot, proci, tag, comm);
                }
            }
        }
        else
        {
            if (x.empty())
            {
                // Nothing to do
            }
            else if (is_contiguous<Type>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    x.cdata_bytes(),
                    x.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream::send(x, UPstream::masterNo(), tag, comm);
            }
        }

        UPstream::waitRequests(startOfRequests);

        // LUBacksubstitute then scatter (globalIndex::scatter pattern)
        if (UPstream::master(comm))
        {
            LUBacksubstitute(coarsestLUPtr_(), coarsestPivots_, allx);

            x = SubList<Type>(allx, x.size());

            for (const int proci : UPstream::subProcs(comm))
            {
                SubList<Type> procSlot
                (
                    allx,
                    coarsestProcOffsets_[proci+1]-coarsestProcOffsets_[proci],
                    coarsestProcOffsets_[proci]
                );

                if (procSlot.empty())
                {
                    // Nothing to do
                }
                else if (is_contiguous<Type>::value)
                {
                    UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        procSlot.cdata_bytes(),
                        procSlot.size_bytes(),
                        tag,
                        comm
                    );
                }
                else
                {
                    OPstream::send(procSlot, proci, tag, comm);
                }
            }
        }
        else
        {
            if (x.empty())
            {
                // Nothing to do
            }
            else if (is_contiguous<Type>::value)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    x.data_bytes(),
                    x.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream::recv(x, UPstream::masterNo(), tag, comm);
            }
        }

        UPstream::waitRequests(startOfRequests);
    }
    else
    {
        LUBacksubstitute(coarsestLUPtr_(), coarsestPivots_, x);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type, class DType, class LUType>
Foam::SolverPerformance<Type>
Foam::LduGAMGSolver<Type, DType, LUType>::solve(Field<Type>& psi) const
{
    SolverPerformance<Type> solverPerf(typeName, this->fieldName_);

    Field<Type> Apsi(psi.size());
    Field<Type> finestCorrection(psi.size());
    Field<Type> finestResidual(psi.size());
    Field<Type> temp(psi.size());

    this->matrix_.Amul(Apsi, psi);

    const Type normFactor = this->normFactor(psi, Apsi, temp);

    finestResidual = this->matrix_.source() - Apsi;

    const auto initialResidual = cmptDivide
    (
        gSumCmptMag(finestResidual),
        normFactor
    );

    if constexpr
    (
        std::is_same<Type, vector4>::value
     && std::is_same<DType, tensor4>::value
    )
    {
        solverPerf.initialResidual() = Type(initialResidual);
    }
    else
    {
        solverPerf.initialResidual() = initialResidual;
    }

    solverPerf.finalResidual() = solverPerf.initialResidual();

    if ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    label nIter = 0;

    if
    (
        this->minIter_ > 0
     || !solverPerf.checkConvergence
        (
            this->tolerance_,
            this->relTol_,
            this->log_
        )
    )
    {
        PtrList<Field<Type>> coarseCorrFields;
        PtrList<Field<Type>> coarseSources;
        PtrList<typename LduMatrix<Type, DType, LUType>::smoother> smoothers;

        initVcycle(coarseCorrFields, coarseSources, smoothers);

        do
        {
            Vcycle
            (
                smoothers,
                psi,
                Apsi,
                finestCorrection,
                finestResidual,
                coarseCorrFields,
                coarseSources
            );

            this->matrix_.Amul(Apsi, psi);
            finestResidual = this->matrix_.source() - Apsi;

            const auto finalResidual = cmptDivide
            (
                gSumCmptMag(finestResidual),
                normFactor
            );

            if constexpr
            (
                std::is_same<Type, vector4>::value
             && std::is_same<DType, tensor4>::value
            )
            {
                solverPerf.finalResidual() = Type(finalResidual);
            }
            else
            {
                solverPerf.finalResidual() = finalResidual;
            }

            if ((this->log_ >= 2) || (LduMatrix<Type, DType, LUType>::debug >= 2))
            {
                solverPerf.print(Info);
            }
        } while
        (
            (
                ++nIter < this->maxIter_
             && !solverPerf.checkConvergence
                (
                    this->tolerance_,
                    this->relTol_,
                    this->log_
                )
            )
         || nIter < this->minIter_
        );
    }

    solverPerf.nIterations() =
        pTraits<typename pTraits<Type>::labelType>::one*nIter;

    return solverPerf;
}


// ************************************************************************* //
