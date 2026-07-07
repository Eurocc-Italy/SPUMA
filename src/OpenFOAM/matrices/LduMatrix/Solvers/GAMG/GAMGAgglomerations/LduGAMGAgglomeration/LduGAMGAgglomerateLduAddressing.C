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

#include "LduGAMGAgglomeration.H"
#include "LduGAMGInterface.H"
#include "processorLduGAMGInterface.H"

void Foam::LduGAMGAgglomeration::agglomerateLduAddressing
(
    const label fineLevelIndex
)
{
    const lduMesh& fineMesh = meshLevel(fineLevelIndex);
    const lduAddressing& fineMeshAddr = fineMesh.lduAddr();

    const labelUList& upperAddr = fineMeshAddr.upperAddr();
    const labelUList& lowerAddr = fineMeshAddr.lowerAddr();

    const label nFineFaces = upperAddr.size();
    const labelField& restrictMap = restrictAddressing(fineLevelIndex);

    if (min(restrictMap) == -1)
    {
        FatalErrorInFunction
            << "min(restrictMap) == -1" << exit(FatalError);
    }

    if (restrictMap.size() != fineMeshAddr.size())
    {
        FatalErrorInFunction
            << "restrict map does not correspond to fine level. " << endl
            << " Sizes: restrictMap: " << restrictMap.size()
            << " nEqns: " << fineMeshAddr.size()
            << abort(FatalError);
    }

    const label nCoarseCells = nCells_[fineLevelIndex];

    label maxNnbrs = 10;
    labelList cCellnFaces(nCoarseCells, Zero);
    labelList cCellFaces(maxNnbrs*nCoarseCells);

    faceRestrictAddressing_.set(fineLevelIndex, new labelList(nFineFaces));
    labelList& faceRestrictAddr = faceRestrictAddressing_[fineLevelIndex];

    labelList initCoarseNeighb(nFineFaces);
    label& nCoarseFaces = nFaces_[fineLevelIndex];
    nCoarseFaces = 0;

    forAll(upperAddr, fineFacei)
    {
        const label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
        const label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

        if (rmUpperAddr == rmLowerAddr)
        {
            faceRestrictAddr[fineFacei] = -(rmUpperAddr + 1);
        }
        else
        {
            label cOwn = rmUpperAddr;
            label cNei = rmLowerAddr;

            if (rmUpperAddr > rmLowerAddr)
            {
                cOwn = rmLowerAddr;
                cNei = rmUpperAddr;
            }

            label* ccFaces = &cCellFaces[maxNnbrs*cOwn];
            bool nbrFound = false;
            label& ccnFaces = cCellnFaces[cOwn];

            for (label i = 0; i < ccnFaces; ++i)
            {
                if (initCoarseNeighb[ccFaces[i]] == cNei)
                {
                    nbrFound = true;
                    faceRestrictAddr[fineFacei] = ccFaces[i];
                    break;
                }
            }

            if (!nbrFound)
            {
                if (ccnFaces >= maxNnbrs)
                {
                    const label oldMaxNnbrs = maxNnbrs;
                    maxNnbrs *= 2;

                    cCellFaces.setSize(maxNnbrs*nCoarseCells);

                    forAllReverse(cCellnFaces, i)
                    {
                        label* oldCcNbrs = &cCellFaces[oldMaxNnbrs*i];
                        label* newCcNbrs = &cCellFaces[maxNnbrs*i];

                        for (label j = 0; j < cCellnFaces[i]; ++j)
                        {
                            newCcNbrs[j] = oldCcNbrs[j];
                        }
                    }

                    ccFaces = &cCellFaces[maxNnbrs*cOwn];
                }

                ccFaces[ccnFaces] = nCoarseFaces;
                initCoarseNeighb[nCoarseFaces] = cNei;
                faceRestrictAddr[fineFacei] = nCoarseFaces;
                ccnFaces++;
                nCoarseFaces++;
            }
        }
    }

    labelList coarseOwner(nCoarseFaces);
    labelList coarseNeighbour(nCoarseFaces);
    labelList coarseFaceMap(nCoarseFaces);

    label coarseFacei = 0;

    forAll(cCellnFaces, coarseCelli)
    {
        label* cFaces = &cCellFaces[maxNnbrs*coarseCelli];
        const label ccnFaces = cCellnFaces[coarseCelli];

        for (label i = 0; i < ccnFaces; ++i)
        {
            coarseOwner[coarseFacei] = coarseCelli;
            coarseNeighbour[coarseFacei] = initCoarseNeighb[cFaces[i]];
            coarseFaceMap[cFaces[i]] = coarseFacei;
            coarseFacei++;
        }
    }

    forAll(faceRestrictAddr, fineFacei)
    {
        if (faceRestrictAddr[fineFacei] >= 0)
        {
            faceRestrictAddr[fineFacei] =
                coarseFaceMap[faceRestrictAddr[fineFacei]];
        }
    }

    faceFlipMap_.set(fineLevelIndex, new boolList(nFineFaces, false));
    boolList& faceFlipMap = faceFlipMap_[fineLevelIndex];

    forAll(faceRestrictAddr, fineFacei)
    {
        const label coarseFacei = faceRestrictAddr[fineFacei];

        if (coarseFacei >= 0)
        {
            const label cOwn = coarseOwner[coarseFacei];
            const label cNei = coarseNeighbour[coarseFacei];

            const label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
            const label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

            if (cOwn == rmUpperAddr && cNei == rmLowerAddr)
            {
                faceFlipMap[fineFacei] = true;
            }
            else if (!(cOwn == rmLowerAddr && cNei == rmUpperAddr))
            {
                FatalErrorInFunction
                    << "problem."
                    << " fineFacei:" << fineFacei
                    << " rmUpperAddr:" << rmUpperAddr
                    << " rmLowerAddr:" << rmLowerAddr
                    << " coarseFacei:" << coarseFacei
                    << " cOwn:" << cOwn
                    << " cNei:" << cNei
                    << exit(FatalError);
            }
        }
    }

    cCellnFaces.clear();
    cCellFaces.clear();
    initCoarseNeighb.clear();
    coarseFaceMap.clear();

    const lduInterfacePtrsList& fineInterfaces = interfaceLevel(fineLevelIndex);

    nPatchFaces_.set
    (
        fineLevelIndex,
        new labelList(fineInterfaces.size(), Zero)
    );
    labelList& nPatchFaces = nPatchFaces_[fineLevelIndex];

    patchFaceRestrictAddressing_.set
    (
        fineLevelIndex,
        new labelListList(fineInterfaces.size())
    );
    labelListList& patchFineToCoarse =
        patchFaceRestrictAddressing_[fineLevelIndex];

    const label startOfRequests = UPstream::nRequests();

    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            if (fineLevelIndex == 0)
            {
                fineInterfaces[inti].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    restrictMap,
                    fineMeshAddr.patchAddr(inti)
                );
            }
            else
            {
                fineInterfaces[inti].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    restrictMap
                );
            }
        }
    }

    UPstream::waitRequests(startOfRequests);

    meshLevels_.set
    (
        fineLevelIndex,
        new lduPrimitiveMesh
        (
            nCoarseCells,
            coarseOwner,
            coarseNeighbour,
            fineMesh.comm(),
            true
        )
    );

    lduInterfacePtrsList coarseInterfaces(fineInterfaces.size());

    forAll(fineInterfaces, inti)
    {
        if (fineInterfaces.set(inti))
        {
            tmp<labelField> restrictMapInternalField;

            if (fineLevelIndex == 0)
            {
                restrictMapInternalField =
                    fineInterfaces[inti].interfaceInternalField
                    (
                        restrictMap,
                        fineMeshAddr.patchAddr(inti)
                    );
            }
            else
            {
                restrictMapInternalField =
                    fineInterfaces[inti].interfaceInternalField(restrictMap);
            }

            tmp<labelField> nbrRestrictMapInternalField =
                fineInterfaces[inti].internalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    restrictMap
                );

            coarseInterfaces.set
            (
                inti,
                LduGAMGInterface::New
                (
                    inti,
                    meshLevels_[fineLevelIndex].rawInterfaces(),
                    fineInterfaces[inti],
                    restrictMapInternalField(),
                    nbrRestrictMapInternalField(),
                    fineLevelIndex,
                    fineMesh.comm()
                ).ptr()
            );

            nPatchFaces[inti] = coarseInterfaces[inti].faceCells().size();
            patchFineToCoarse[inti] = refCast<const LduGAMGInterface>
            (
                coarseInterfaces[inti]
            ).faceRestrictAddressing();
        }
    }

    meshLevels_[fineLevelIndex].addInterfaces
    (
        coarseInterfaces,
        lduPrimitiveMesh::nonBlockingSchedule<processorLduGAMGInterface>
        (
            coarseInterfaces
        )
    );
}


void Foam::LduGAMGAgglomeration::combineLevels(const label curLevel)
{
    const label prevLevel = curLevel - 1;

    nCells_[prevLevel] = nCells_[curLevel];
    nFaces_[prevLevel] = nFaces_[curLevel];

    const labelList& curResAddr = restrictAddressing_[curLevel];
    labelList& prevResAddr = restrictAddressing_[prevLevel];

    const labelList& curFaceResAddr = faceRestrictAddressing_[curLevel];
    labelList& prevFaceResAddr = faceRestrictAddressing_[prevLevel];
    const boolList& curFaceFlipMap = faceFlipMap_[curLevel];
    boolList& prevFaceFlipMap = faceFlipMap_[prevLevel];

    forAll(prevFaceResAddr, i)
    {
        if (prevFaceResAddr[i] >= 0)
        {
            const label fineFacei = prevFaceResAddr[i];
            prevFaceResAddr[i] = curFaceResAddr[fineFacei];
            prevFaceFlipMap[i] = curFaceFlipMap[fineFacei];
        }
        else
        {
            const label fineFacei = -prevFaceResAddr[i] - 1;
            prevFaceResAddr[i] = -curResAddr[fineFacei] - 1;
            prevFaceFlipMap[i] = curFaceFlipMap[fineFacei];
        }
    }

    faceRestrictAddressing_.set(curLevel, nullptr);
    faceFlipMap_.set(curLevel, nullptr);

    forAll(prevResAddr, i)
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    const labelListList& curPatchFaceResAddr =
        patchFaceRestrictAddressing_[curLevel];
    labelListList& prevPatchFaceResAddr =
        patchFaceRestrictAddressing_[prevLevel];

    forAll(prevPatchFaceResAddr, inti)
    {
        const labelList& curPatchResAddr = curPatchFaceResAddr[inti];
        labelList& prevPatchResAddr = prevPatchFaceResAddr[inti];

        forAll(prevPatchResAddr, i)
        {
            const label fineFacei = prevPatchResAddr[i];
            prevPatchResAddr[i] = curPatchResAddr[fineFacei];
        }
    }

    restrictAddressing_.set(curLevel, nullptr);
    nPatchFaces_[prevLevel] = nPatchFaces_[curLevel];

    const lduInterfacePtrsList& curInterLevel =
        meshLevels_[curLevel].rawInterfaces();
    const lduInterfacePtrsList& prevInterLevel =
        meshLevels_[prevLevel].rawInterfaces();

    forAll(prevInterLevel, inti)
    {
        if (prevInterLevel.set(inti))
        {
            LduGAMGInterface& prevInt = refCast<LduGAMGInterface>
            (
                const_cast<lduInterface&>(prevInterLevel[inti])
            );
            const LduGAMGInterface& curInt = refCast<const LduGAMGInterface>
            (
                curInterLevel[inti]
            );

            prevInt.combine(curInt);
        }
    }

    meshLevels_.set(prevLevel, meshLevels_.set(curLevel, nullptr));
}


// ************************************************************************* //