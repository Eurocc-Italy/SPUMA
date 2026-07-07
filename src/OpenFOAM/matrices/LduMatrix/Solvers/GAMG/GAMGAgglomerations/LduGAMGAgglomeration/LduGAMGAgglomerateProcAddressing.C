/*---------------------------------------------------------------------------*\
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

Description
    Native processor agglomeration for the templated LduMatrix GAMG path.
    Ports the gather/combine addressing logic from the scalar GAMG and
    lduPrimitiveMesh, producing native LduGAMGInterface geometry.

\*---------------------------------------------------------------------------*/

#include "LduGAMGAgglomeration.H"
#include "LduGAMGInterface.H"
#include "processorLduGAMGInterface.H"
#include "processorLduInterface.H"
#include "PstreamBuffers.H"
#include "globalIndex.H"
#include "CompactListList.H"
#include "PrecisionAdaptor.H"
#include "edgeHashes.H"
#include "Map.H"
#include "ListOps.H"

// * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * * //

void Foam::LduGAMGAgglomeration::calculateRegionMaster
(
    const label comm,
    const labelList& procAgglomMap,
    labelList& masterProcs,
    List<label>& agglomProcIDs
)
{
    // Determine the master processors
    Map<label> agglomToMaster(procAgglomMap.size());

    forAll(procAgglomMap, proci)
    {
        const label coarsei = procAgglomMap[proci];

        auto iter = agglomToMaster.find(coarsei);
        if (iter.good())
        {
            iter.val() = min(iter.val(), proci);
        }
        else
        {
            agglomToMaster.insert(coarsei, proci);
        }
    }

    masterProcs.setSize(agglomToMaster.size());
    forAllConstIters(agglomToMaster, iter)
    {
        masterProcs[iter.key()] = iter.val();
    }


    // Collect all the processors in my agglomeration
    const label myProcID = UPstream::myProcNo(comm);
    const label myAgglom = procAgglomMap[myProcID];

    // Get all processors agglomerating to the same coarse processor
    agglomProcIDs = findIndices(procAgglomMap, myAgglom);

    // Make sure the master is the first element.
    const label index = agglomProcIDs.find(agglomToMaster[myAgglom]);

    std::swap(agglomProcIDs[0], agglomProcIDs[index]);
}


Foam::label Foam::LduGAMGAgglomeration::findConnectedInterface
(
    const lduMesh& myMesh,
    const PtrList<lduPrimitiveMesh>& otherMeshes,
    const labelPairList& procAndInterfaces,
    const label nbrProci,
    const label myRank
)
{
    label nbrInti = -1;
    for (const auto& procAndInterface : procAndInterfaces)
    {
        const label proci = procAndInterface[0];

        if (proci == nbrProci)
        {
            const label interfacei = procAndInterface[1];
            const lduInterfacePtrsList interfaces =
                lduPrimitiveMesh::mesh(myMesh, otherMeshes, proci).interfaces();
            const processorLduInterface& pldui =
                refCast<const processorLduInterface>(interfaces[interfacei]);

            if (pldui.neighbProcNo() == myRank)
            {
                nbrInti = procAndInterface[1];
                break;
            }
        }
    }

    if (nbrInti == -1)
    {
        FatalErrorInFunction
            << "procAndInterfaces:" << procAndInterfaces << abort(FatalError);
    }
    return nbrInti;
}


void Foam::LduGAMGAgglomeration::gatherMeshes
(
    const label comm,
    const lduMesh& mesh,
    PtrList<lduPrimitiveMesh>& otherMeshes
)
{
    // Force calculation of schedule (since does parallel comms)
    (void)mesh.lduAddr().patchSchedule();

    PstreamBuffers pBufs(comm);

    // Send to master
    if (!UPstream::master(comm))
    {
        const lduAddressing& addressing = mesh.lduAddr();

        lduInterfacePtrsList interfaces(mesh.interfaces());
        boolList validInterface(interfaces.size());
        forAll(interfaces, intI)
        {
            validInterface[intI] = interfaces.set(intI);
        }

        UOPstream toMaster(UPstream::masterNo(), pBufs);

        toMaster
            << addressing.size()
            << addressing.lowerAddr()
            << addressing.upperAddr()
            << validInterface;

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const LduGAMGInterface& interface =
                    refCast<const LduGAMGInterface>(interfaces[intI]);

                toMaster << interface.type();
                interface.write(toMaster);
            }
        }
    }

    pBufs.finishedGathers();

    // Consume
    if (UPstream::master(comm))
    {
        const label nProcs = UPstream::nProcs(comm);

        otherMeshes.setSize(nProcs - 1);

        for (const int proci : UPstream::subProcs(comm))
        {
            UIPstream fromProc(proci, pBufs);

            const label nCells = readLabel(fromProc);
            labelList lowerAddr(fromProc);
            labelList upperAddr(fromProc);
            const boolList validInterface(fromProc);

            // Construct mesh without interfaces
            otherMeshes.set
            (
                proci - 1,
                new lduPrimitiveMesh
                (
                    nCells,
                    lowerAddr,
                    upperAddr,
                    mesh.comm(),
                    true    // reuse
                )
            );

            // Construct native GAMG interfaces
            lduInterfacePtrsList newInterfaces(validInterface.size());
            forAll(validInterface, intI)
            {
                if (validInterface[intI])
                {
                    word coupleType(fromProc);

                    newInterfaces.set
                    (
                        intI,
                        LduGAMGInterface::New
                        (
                            coupleType,
                            intI,
                            otherMeshes[proci - 1].rawInterfaces(),
                            fromProc
                        ).ptr()
                    );
                }
            }

            otherMeshes[proci - 1].addInterfaces
            (
                newInterfaces,
                lduPrimitiveMesh::nonBlockingSchedule<processorLduGAMGInterface>
                (
                    newInterfaces
                )
            );
        }
    }
}


Foam::autoPtr<Foam::lduPrimitiveMesh>
Foam::LduGAMGAgglomeration::combineMeshes
(
    const label comm,
    const labelList& procAgglomMap,
    const labelList& procIDs,
    const lduMesh& myMesh,
    const PtrList<lduPrimitiveMesh>& otherMeshes,
    labelList& cellOffsets,
    labelList& faceOffsets,
    labelListList& faceMap,
    labelListList& boundaryMap,
    labelListListList& boundaryFaceMap
)
{
    const label currentComm = myMesh.comm();
    const label nMeshes = otherMeshes.size() + 1;

    const label myAgglom = procAgglomMap[UPstream::myProcNo(currentComm)];

    const lduInterfacePtrsList myInterfaces = myMesh.interfaces();

    forAll(procIDs, i)
    {
        if (procAgglomMap[procIDs[i]] != procAgglomMap[procIDs[0]])
        {
            FatalErrorInFunction
                << "Processor " << procIDs[i]
                << " agglomerates to " << procAgglomMap[procIDs[i]]
                << " whereas other processors " << procIDs
                << " agglomerate to "
                << labelUIndList(procAgglomMap, procIDs)
                << exit(FatalError);
        }
    }


    // Cells get added in order.
    cellOffsets.setSize(nMeshes + 1);
    cellOffsets[0] = 0;
    for (label procMeshI = 0; procMeshI < nMeshes; ++procMeshI)
    {
        const lduMesh& procMesh =
            lduPrimitiveMesh::mesh(myMesh, otherMeshes, procMeshI);

        cellOffsets[procMeshI + 1] =
            cellOffsets[procMeshI] + procMesh.lduAddr().size();
    }


    // Merged proc interfaces: map from two coarse processors back to
    // (procMesh, interface)
    EdgeMap<labelPairList> mergedMap(2*myInterfaces.size());

    // Unmerged (surviving) proc interfaces
    EdgeMap<labelPairList> unmergedMap(mergedMap.size());

    // (unmerged) global/non-processor interfaces - unsupported here
    List<DynamicList<label>> procToGlobal(nMeshes);


    boundaryMap.setSize(nMeshes);
    boundaryFaceMap.setSize(nMeshes);

    labelList nCoupledFaces(nMeshes, Zero);

    for (label procMeshI = 0; procMeshI < nMeshes; ++procMeshI)
    {
        const lduInterfacePtrsList interfaces =
            lduPrimitiveMesh::mesh(myMesh, otherMeshes, procMeshI).interfaces();

        boundaryMap[procMeshI].setSize(interfaces.size(), -1);
        boundaryFaceMap[procMeshI].setSize(interfaces.size());

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const lduInterface& ldui = interfaces[intI];

                if (isA<processorLduInterface>(ldui))
                {
                    const processorLduInterface& pldui =
                        refCast<const processorLduInterface>(ldui);

                    const label agglom0 = procAgglomMap[pldui.myProcNo()];
                    const label agglom1 = procAgglomMap[pldui.neighbProcNo()];

                    const edge procEdge(agglom0, agglom1);

                    if (agglom0 != myAgglom && agglom1 != myAgglom)
                    {
                        FatalErrorInFunction
                            << "At mesh from processor " << procIDs[procMeshI]
                            << " have interface " << intI
                            << " with myProcNo:" << pldui.myProcNo()
                            << " with neighbProcNo:" << pldui.neighbProcNo()
                            << exit(FatalError);
                    }
                    else if (agglom0 == myAgglom && agglom1 == myAgglom)
                    {
                        // Merged interface
                        const label nbrProcMeshI =
                            procIDs.find(pldui.neighbProcNo());

                        if (procMeshI < nbrProcMeshI)
                        {
                            nCoupledFaces[procMeshI] += ldui.faceCells().size();
                        }

                        mergedMap(procEdge).append
                        (
                            labelPair(procMeshI, intI)
                        );
                    }
                    else
                    {
                        unmergedMap(procEdge).append
                        (
                            labelPair(procMeshI, intI)
                        );
                    }
                }
                else
                {
                    // Non-processor (global/cyclic) interface - unsupported
                    procToGlobal[procMeshI].append(intI);
                }
            }
        }
    }


    if (procToGlobal[0].size())
    {
        FatalErrorInFunction
            << "Non-processor coupled interfaces surviving at the"
            << " processor-agglomerated level are not supported in the"
            << " native LduMatrix GAMG path." << nl
            << "Number of such interfaces on master mesh: "
            << procToGlobal[0].size()
            << exit(FatalError);
    }


    // Adapt faceOffsets for internal interfaces
    faceOffsets.setSize(nMeshes + 1);
    faceOffsets[0] = 0;
    faceMap.setSize(nMeshes);
    for (label procMeshI = 0; procMeshI < nMeshes; ++procMeshI)
    {
        const lduMesh& procMesh =
            lduPrimitiveMesh::mesh(myMesh, otherMeshes, procMeshI);
        const label nInternal = procMesh.lduAddr().lowerAddr().size();

        faceOffsets[procMeshI + 1] =
            faceOffsets[procMeshI] + nInternal + nCoupledFaces[procMeshI];

        labelList& map = faceMap[procMeshI];
        map.setSize(nInternal);
        forAll(map, i)
        {
            map[i] = faceOffsets[procMeshI] + i;
        }
    }


    // Combine upper and lower
    labelList lowerAddr(faceOffsets.last(), -1);
    labelList upperAddr(lowerAddr.size(), -1);


    // Old internal faces and resolved coupled interfaces
    for (label procMeshI = 0; procMeshI < nMeshes; ++procMeshI)
    {
        const lduMesh& procMesh =
            lduPrimitiveMesh::mesh(myMesh, otherMeshes, procMeshI);

        const labelUList& l = procMesh.lduAddr().lowerAddr();
        const labelUList& u = procMesh.lduAddr().upperAddr();

        label allFacei = faceOffsets[procMeshI];

        forAll(l, facei)
        {
            lowerAddr[allFacei] = cellOffsets[procMeshI] + l[facei];
            upperAddr[allFacei] = cellOffsets[procMeshI] + u[facei];
            allFacei++;
        }


        const lduInterfacePtrsList interfaces = procMesh.interfaces();

        forAll(interfaces, intI)
        {
            if (interfaces.set(intI))
            {
                const lduInterface& ldui = interfaces[intI];

                if (isA<processorLduInterface>(ldui))
                {
                    const processorLduInterface& pldui =
                        refCast<const processorLduInterface>(ldui);

                    const label myP = pldui.myProcNo();
                    const label nbrP = pldui.neighbProcNo();
                    const label nbrProcMeshI = procIDs.find(nbrP);

                    if (procMeshI < nbrProcMeshI)
                    {
                        const label agglom0 = procAgglomMap[myP];
                        const label agglom1 = procAgglomMap[nbrP];

                        const auto fnd =
                            mergedMap.cfind(edge(agglom0, agglom1));

                        if (fnd.good())
                        {
                            const labelPairList& elems = fnd();

                            const label nbrIntI = findConnectedInterface
                            (
                                myMesh,
                                otherMeshes,
                                elems,
                                nbrProcMeshI,
                                procIDs[procMeshI]
                            );

                            const lduInterfacePtrsList nbrInterfaces =
                                lduPrimitiveMesh::mesh
                                (
                                    myMesh,
                                    otherMeshes,
                                    nbrProcMeshI
                                ).interfaces();

                            const labelUList& faceCells = ldui.faceCells();
                            const labelUList& nbrFaceCells =
                                nbrInterfaces[nbrIntI].faceCells();

                            if (faceCells.size() != nbrFaceCells.size())
                            {
                                FatalErrorInFunction
                                    << "faceCells:" << faceCells
                                    << " nbrFaceCells:" << nbrFaceCells
                                    << abort(FatalError);
                            }

                            labelList& bfMap =
                                boundaryFaceMap[procMeshI][intI];
                            labelList& nbrBfMap =
                                boundaryFaceMap[nbrProcMeshI][nbrIntI];

                            bfMap.setSize(faceCells.size());
                            nbrBfMap.setSize(faceCells.size());

                            forAll(faceCells, pfI)
                            {
                                lowerAddr[allFacei] =
                                    cellOffsets[procMeshI] + faceCells[pfI];
                                bfMap[pfI] = allFacei;
                                upperAddr[allFacei] =
                                    cellOffsets[nbrProcMeshI]
                                  + nbrFaceCells[pfI];
                                nbrBfMap[pfI] = (-allFacei - 1);
                                allFacei++;
                            }
                        }
                    }
                }
            }
        }
    }


    // Sort upper-tri order
    {
        labelList oldToNew
        (
            lduPrimitiveMesh::upperTriOrder
            (
                cellOffsets.last(),
                lowerAddr,
                upperAddr
            )
        );

        forAll(faceMap, procMeshI)
        {
            labelList& map = faceMap[procMeshI];
            forAll(map, i)
            {
                if (map[i] >= 0)
                {
                    map[i] = oldToNew[map[i]];
                }
                else
                {
                    const label allFacei = -map[i] - 1;
                    map[i] = -oldToNew[allFacei] - 1;
                }
            }
        }

        inplaceReorder(oldToNew, lowerAddr);
        inplaceReorder(oldToNew, upperAddr);

        forAll(boundaryFaceMap, proci)
        {
            const labelList& bMap = boundaryMap[proci];
            forAll(bMap, intI)
            {
                if (bMap[intI] == -1)
                {
                    labelList& bfMap = boundaryFaceMap[proci][intI];

                    forAll(bfMap, i)
                    {
                        if (bfMap[i] >= 0)
                        {
                            bfMap[i] = oldToNew[bfMap[i]];
                        }
                        else
                        {
                            const label allFacei = -bfMap[i] - 1;
                            bfMap[i] = (-oldToNew[allFacei] - 1);
                        }
                    }
                }
            }
        }
    }


    // Construct the combined mesh (without interfaces yet)
    const label nCells = cellOffsets.last();

    autoPtr<lduPrimitiveMesh> meshPtr
    (
        new lduPrimitiveMesh
        (
            nCells,
            lowerAddr,
            upperAddr,
            comm,
            true    // reuse
        )
    );
    lduPrimitiveMesh& newMesh = meshPtr();


    // Kept (unmerged) processor interfaces
    lduInterfacePtrsList newInterfaces(unmergedMap.size());

    label allInterfacei = 0;

    for (const auto& iter : unmergedMap.csorted())
    {
        const labelPairList& elems = iter.val();

        // Sort processors in increasing order so both sides walk in same order
        labelPairList procPairs(elems.size());
        forAll(elems, i)
        {
            const labelPair& elem = elems[i];
            const label procMeshI = elem[0];
            const label interfacei = elem[1];
            const lduInterfacePtrsList interfaces =
                lduPrimitiveMesh::mesh
                (
                    myMesh,
                    otherMeshes,
                    procMeshI
                ).interfaces();

            const processorLduInterface& pldui =
                refCast<const processorLduInterface>(interfaces[interfacei]);
            const label myProcNo = pldui.myProcNo();
            const label nbrProcNo = pldui.neighbProcNo();
            procPairs[i] = labelPair
            (
                min(myProcNo, nbrProcNo),
                max(myProcNo, nbrProcNo)
            );
        }

        labelList order(sortedOrder(procPairs));

        // Count
        label n = 0;
        forAll(order, i)
        {
            const labelPair& elem = elems[order[i]];
            const label procMeshI = elem[0];
            const label interfacei = elem[1];
            const lduInterfacePtrsList interfaces =
                lduPrimitiveMesh::mesh
                (
                    myMesh,
                    otherMeshes,
                    procMeshI
                ).interfaces();

            n += interfaces[interfacei].faceCells().size();
        }

        // Size and fill
        labelField allFaceCells(n);
        labelField allFaceRestrictAddressing(n);
        n = 0;

        forAll(order, i)
        {
            const labelPair& elem = elems[order[i]];
            const label procMeshI = elem[0];
            const label interfacei = elem[1];
            const lduInterfacePtrsList interfaces =
                lduPrimitiveMesh::mesh
                (
                    myMesh,
                    otherMeshes,
                    procMeshI
                ).interfaces();

            boundaryMap[procMeshI][interfacei] = allInterfacei;
            labelList& bfMap = boundaryFaceMap[procMeshI][interfacei];

            const labelUList& l = interfaces[interfacei].faceCells();
            bfMap.setSize(l.size());

            forAll(l, facei)
            {
                allFaceCells[n] = cellOffsets[procMeshI] + l[facei];
                allFaceRestrictAddressing[n] = n;
                bfMap[facei] = n;
                n++;
            }
        }


        // Find out remote processor in new communicator
        label neighbProcNo = -1;

        if (iter.key()[0] == myAgglom)
        {
            if (iter.key()[1] == myAgglom)
            {
                FatalErrorInFunction
                    << "problem procEdge:" << iter.key()
                    << exit(FatalError);
            }
            neighbProcNo = iter.key()[1];
        }
        else
        {
            if (iter.key()[1] != myAgglom)
            {
                FatalErrorInFunction
                    << "problem procEdge:" << iter.key()
                    << exit(FatalError);
            }
            neighbProcNo = iter.key()[0];
        }

        newInterfaces.set
        (
            allInterfacei,
            new processorLduGAMGInterface
            (
                allInterfacei,
                newMesh.rawInterfaces(),
                allFaceCells,
                allFaceRestrictAddressing,
                comm,
                myAgglom,
                neighbProcNo,
                tensorField(),          // forwardT
                UPstream::msgType()     // tag
            )
        );

        allInterfacei++;
    }


    newMesh.addInterfaces
    (
        newInterfaces,
        lduPrimitiveMesh::nonBlockingSchedule<processorLduGAMGInterface>
        (
            newInterfaces
        )
    );

    return meshPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LduGAMGAgglomeration::clearLevel(const label i)
{
    if (hasMeshLevel(i))
    {
        meshLevels_.set(i - 1, nullptr);

        if (i < nCells_.size())
        {
            nCells_[i] = -555;
            restrictAddressing_.set(i, nullptr);
            nFaces_[i] = -666;
            faceRestrictAddressing_.set(i, nullptr);
            faceFlipMap_.set(i, nullptr);
            nPatchFaces_.set(i, nullptr);
            patchFaceRestrictAddressing_.set(i, nullptr);
        }
    }
}


void Foam::LduGAMGAgglomeration::procAgglomerateLduAddressing
(
    const label meshComm,
    const labelList& procAgglomMap,
    const labelList& procIDs,
    const label allMeshComm,
    const label levelIndex
)
{
    const lduMesh& myMesh = meshLevels_[levelIndex - 1];
    const label nOldInterfaces = myMesh.interfaces().size();

    procAgglomMap_.set(levelIndex, new labelList(procAgglomMap));
    agglomProcIDs_.set(levelIndex, new labelList(procIDs));
    procCommunicator_[levelIndex] = allMeshComm;

    procAgglomCommunicator_.set
    (
        levelIndex,
        new UPstream::communicator(meshComm, procIDs)
    );
    const label comm = agglomCommunicator(levelIndex);

    // Default-initialise (also kept on slaves for convenience)
    procCellOffsets_.set(levelIndex, new labelList(0));
    procFaceMap_.set(levelIndex, new labelListList(0));
    procBoundaryMap_.set(levelIndex, new labelListList(0));
    procBoundaryFaceMap_.set(levelIndex, new labelListListList(0));


    // Collect meshes onto master
    PtrList<lduPrimitiveMesh> otherMeshes;
    gatherMeshes(comm, myMesh, otherMeshes);

    if (UPstream::myProcNo(meshComm) == procIDs[0])
    {
        labelList procFaceOffsets;

        meshLevels_.set
        (
            levelIndex - 1,
            combineMeshes
            (
                allMeshComm,
                procAgglomMap,
                procIDs,
                myMesh,
                otherMeshes,

                procCellOffsets_[levelIndex],
                procFaceOffsets,
                procFaceMap_[levelIndex],
                procBoundaryMap_[levelIndex],
                procBoundaryFaceMap_[levelIndex]
            ).ptr()
        );
    }


    // Scatter the procBoundaryMap back to the originating processor so it
    // knows which proc boundaries are kept (used to limit interfaceField
    // transfer in procAgglomerateMatrix).
    {
        const CompactListList<label> data
        (
            CompactListList<label>::pack<labelList>
            (
                procBoundaryMap_[levelIndex]
            )
        );
        const labelList localSizes = data.localSizes();
        const labelList& localStarts = data.offsets();

        procBoundaryMap_[levelIndex].setSize(procIDs.size());
        labelList& bMap =
            procBoundaryMap_[levelIndex][UPstream::myProcNo(comm)];
        bMap.setSize(nOldInterfaces);

        UPstream::scatter
        (
            data.values().cdata(),

            ConstPrecisionAdaptor<int, label, List>(localSizes).cref(),
            ConstPrecisionAdaptor<int, label, List>(localStarts).cref(),

            bMap.data(),
            bMap.size(),
            comm
        );
    }


    // Combine restrict addressing
    procAgglomerateRestrictAddressing(meshComm, procIDs, levelIndex);

    if (UPstream::myProcNo(meshComm) != procIDs[0])
    {
        clearLevel(levelIndex);
    }
}


void Foam::LduGAMGAgglomeration::procAgglomerateRestrictAddressing
(
    const label comm,
    const labelList& procIDs,
    const label levelIndex
)
{
    const bool master =
    (
        UPstream::myProcNo(comm) == (procIDs.empty() ? 0 : procIDs[0])
    );

    labelList fineOffsets;
    labelList coarseOffsets;

    List<labelPair> sizes = globalIndex::listGatherValues
    (
        comm,
        procIDs,
        labelPair
        (
            restrictAddressing_[levelIndex].size(),
            nCells_[levelIndex]
        ),
        UPstream::msgType(),
        UPstream::commsTypes::scheduled
    );

    if (master && !sizes.empty())
    {
        const label len = sizes.size();

        fineOffsets.resize(len+1);
        coarseOffsets.resize(len+1);

        label fineCount = 0;
        label coarseCount = 0;

        for (label i = 0; i < len; ++i)
        {
            fineOffsets[i] = fineCount;
            fineCount += sizes[i].first();

            coarseOffsets[i] = coarseCount;
            coarseCount += sizes[i].second();
        }

        fineOffsets[len] = fineCount;
        coarseOffsets[len] = coarseCount;
    }

    // (cell)restrictAddressing
    labelList procRestrictAddressing;
    if (master)
    {
        procRestrictAddressing.resize(fineOffsets.back());
    }

    globalIndex::gather
    (
        fineOffsets,
        comm,
        procIDs,
        restrictAddressing_[levelIndex],
        procRestrictAddressing,
        UPstream::msgType(),
        UPstream::commsTypes::nonBlocking
    );


    if (master)
    {
        nCells_[levelIndex] = coarseOffsets.last();

        // Renumber consecutively
        for (label proci = 1; proci < procIDs.size(); ++proci)
        {
            SubList<label> procSlot
            (
                procRestrictAddressing,
                fineOffsets[proci+1]-fineOffsets[proci],
                fineOffsets[proci]
            );

            forAll(procSlot, i)
            {
                procSlot[i] += coarseOffsets[proci];
            }
        }

        restrictAddressing_[levelIndex].transfer(procRestrictAddressing);
    }
}


void Foam::LduGAMGAgglomeration::procAgglomerate()
{
    if (size() < 1)
    {
        return;
    }

    // Agglomerate one but last level (also agglomerating restrictAddressing)
    const label fineLevelIndex = size() - 1;

    if (fineLevelIndex <= 0 || !hasMeshLevel(fineLevelIndex))
    {
        return;
    }

    const lduMesh& levelMesh = meshLevel(fineLevelIndex);
    const label levelComm = levelMesh.comm();
    const label nProcs = UPstream::nProcs(levelComm);

    if (nProcs <= 1)
    {
        return;
    }

    // Processor restriction map: per processor the coarse processor
    labelList procAgglomMap(nProcs);
    if (nProcessorsPerMaster_ > 0)
    {
        forAll(procAgglomMap, fineProci)
        {
            procAgglomMap[fineProci] = fineProci / nProcessorsPerMaster_;
        }
    }
    else
    {
        procAgglomMap = Zero;
    }

    labelList masterProcs;
    List<label> agglomProcIDs;
    calculateRegionMaster(levelComm, procAgglomMap, masterProcs, agglomProcIDs);

    // Communicator for the processor-agglomerated matrix (the masters)
    allocatedComms_.push_back
    (
        UPstream::communicator(levelComm, masterProcs)
    );
    const label procAgglomComm = allocatedComms_.back().comm();

    if (UPstream::myProcNo(levelComm) != -1)
    {
        // Collect meshes and restrictAddressing onto master. Overwrites the
        // fine mesh (meshLevels_[fineLevelIndex-1]) and the addressing from
        // fine mesh to coarse mesh (restrictAddressing_[fineLevelIndex]).
        procAgglomerateLduAddressing
        (
            levelComm,
            procAgglomMap,
            agglomProcIDs,
            procAgglomComm,
            fineLevelIndex
        );

        // Combine deeper restrict addressing onto master
        for
        (
            label levelI = fineLevelIndex + 1;
            levelI < meshLevels_.size();
            levelI++
        )
        {
            procAgglomerateRestrictAddressing
            (
                levelComm,
                agglomProcIDs,
                levelI
            );
        }

        if (UPstream::myProcNo(levelComm) == agglomProcIDs[0])
        {
            // On master: recreate coarse meshes from restrict addressing
            for
            (
                label levelI = fineLevelIndex;
                levelI < meshLevels_.size();
                levelI++
            )
            {
                agglomerateLduAddressing(levelI);
            }
        }
        else
        {
            // Agglomerated away: clear mesh storage
            for
            (
                label levelI = fineLevelIndex + 1;
                levelI <= size();
                levelI++
            )
            {
                clearLevel(levelI);
            }
        }
    }
}


// ************************************************************************* //
