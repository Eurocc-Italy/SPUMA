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

Class
    Foam::LduGAMGAgglomeration

Description
    Native serial GAMG agglomeration hierarchy for coupled LduMatrix solvers.

SourceFiles
    LduGAMGAgglomeration.C
    LduGAMGAgglomerationTemplates.C
    LduGAMGAgglomerateLduAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LduGAMGAgglomeration_H
#define Foam_LduGAMGAgglomeration_H

#include "LduMatrix.H"
#include "boolList.H"
#include "lduInterfacePtrsList.H"
#include "lduPrimitiveMesh.H"
#include "primitiveFields.H"
#include "UPstream.H"
#include "DynamicList.H"

namespace Foam
{
class LduGAMGAgglomeration
{
protected:

    // Protected Data

        //- Maximum number of agglomerated levels.
        const label maxLevels_;

        //- Number of cells in the coarsest level.
        label nCellsInCoarsestLevel_;

        //- Finest-level mesh.
        const lduMesh& mesh_;

        //- Finest-level interfaces.
        const lduInterfacePtrsList meshInterfaces_;

        //- Number of coarse cells for each fine level.
        labelList nCells_;

        //- Fine-to-coarse cell restriction addressing.
        PtrList<labelField> restrictAddressing_;

        //- Number of coarse faces for each fine level.
        labelList nFaces_;

        //- Fine-to-coarse face restriction addressing.
        PtrList<labelList> faceRestrictAddressing_;

        //- Face orientation map for restricted faces.
        PtrList<boolList> faceFlipMap_;

        //- Number of coarse patch faces for each fine level.
        PtrList<labelList> nPatchFaces_;

        //- Fine-to-coarse patch-face restriction addressing.
        PtrList<labelListList> patchFaceRestrictAddressing_;

        //- Coarse mesh hierarchy, indexed by fine level.
        PtrList<lduPrimitiveMesh> meshLevels_;


    // Processor agglomeration data (parallel only)

        //- Enable processor agglomeration of the coarsest level.
        bool doProcessorAgglomerate_;

        //- Number of fine processors agglomerated onto each master.
        //  0 means agglomerate everything onto a single master.
        label nProcessorsPerMaster_;

        //- For each level: fine-to-coarse processor restriction.
        mutable PtrList<labelList> procAgglomMap_;

        //- For each level: the processors agglomerated onto the master.
        mutable PtrList<labelList> agglomProcIDs_;

        //- For each level: the agglomerated mesh communicator (-1 if none).
        mutable labelList procCommunicator_;

        //- For each level: communicator spanning the agglomerated processors.
        mutable PtrList<UPstream::communicator> procAgglomCommunicator_;

        //- For each level: cell offsets per agglomerated processor.
        mutable PtrList<labelList> procCellOffsets_;

        //- For each level: internal-face map per agglomerated processor.
        mutable PtrList<labelListList> procFaceMap_;

        //- For each level: boundary map per agglomerated processor.
        mutable PtrList<labelListList> procBoundaryMap_;

        //- For each level: boundary-face map per agglomerated processor.
        mutable PtrList<labelListListList> procBoundaryFaceMap_;

        //- Persistent storage for allocated master communicators.
        mutable DynamicList<UPstream::communicator> allocatedComms_;


    // Protected Member Functions

        //- Construct from finest mesh and controls.
        LduGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Return the selected native agglomerator for the given face weights.
        static autoPtr<LduGAMGAgglomeration> New
        (
            const lduMesh& mesh,
            const scalarField& faceWeights,
            const dictionary& controlDict
        );

        //- Assemble coarse mesh addressing.
        void agglomerateLduAddressing(const label fineLevelIndex);

        //- Combine a level with the previous one.
        void combineLevels(const label curLevel);

        //- Shrink the hierarchy storage to the created levels.
        void compactLevels
        (
            const label nCreatedLevels,
            const bool doProcessorAgglomerate = false
        );

        //- Clear a single level (used on non-master after agglomeration).
        void clearLevel(const label leveli);

        //- Drive processor agglomeration of the coarsest level (masterCoarsest).
        void procAgglomerate();

        //- Gather meshes from a group of processors onto the master.
        static void gatherMeshes
        (
            const label comm,
            const lduMesh& mesh,
            PtrList<lduPrimitiveMesh>& otherMeshes
        );

        //- Combine gathered meshes into a single mesh on the master.
        static autoPtr<lduPrimitiveMesh> combineMeshes
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
        );

        //- Find the interface on nbrProci connected to myRank.
        static label findConnectedInterface
        (
            const lduMesh& myMesh,
            const PtrList<lduPrimitiveMesh>& otherMeshes,
            const labelPairList& procAndInterfaces,
            const label nbrProci,
            const label myRank
        );

        //- Gather and combine the LDU addressing of a level onto the master.
        void procAgglomerateLduAddressing
        (
            const label meshComm,
            const labelList& procAgglomMap,
            const labelList& procIDs,
            const label allMeshComm,
            const label levelIndex
        );

        //- Gather and renumber restrict addressing of a level onto the master.
        void procAgglomerateRestrictAddressing
        (
            const label comm,
            const labelList& procIDs,
            const label levelIndex
        );

        //- Determine the master processors for a processor restriction.
        static void calculateRegionMaster
        (
            const label comm,
            const labelList& procAgglomMap,
            labelList& masterProcs,
            List<label>& agglomProcIDs
        );

        //- Check if further agglomeration should continue.
        bool continueAgglomerating
        (
            const label nCellsInCoarsestLevel,
            const label nCells,
            const label nCoarseCells,
            const label comm
        ) const;

        //- Coarsest-level cell target.
        label nCellsInCoarsestLevel() const
        {
            return nCellsInCoarsestLevel_;
        }

public:

    virtual ~LduGAMGAgglomeration() = default;

    template<class Type, class DType, class LUType>
    static autoPtr<LduGAMGAgglomeration> New
    (
        const LduMatrix<Type, DType, LUType>& matrix,
        const dictionary& controlDict
    );

    label size() const
    {
        return meshLevels_.size();
    }

    //- Return the LDU mesh at a given level.
    const lduMesh& meshLevel(const label leveli) const;

    //- Return true if the level exists.
    bool hasMeshLevel(const label leveli) const;

    //- Return the interface addressing at a given level.
    const lduInterfacePtrsList& interfaceLevel(const label leveli) const;

    //- Return cell restriction addressing for a fine level.
    const labelField& restrictAddressing(const label leveli) const
    {
        return restrictAddressing_[leveli];
    }

    //- Return face restriction addressing for a fine level.
    const labelList& faceRestrictAddressing(const label leveli) const
    {
        return faceRestrictAddressing_[leveli];
    }

    //- Return patch-face restriction addressing for a fine level.
    const labelListList& patchFaceRestrictAddressing(const label leveli) const
    {
        return patchFaceRestrictAddressing_[leveli];
    }

    //- Return the face flip map for a fine level.
    const boolList& faceFlipMap(const label leveli) const
    {
        return faceFlipMap_[leveli];
    }

    //- Return the number of patch faces for a fine level.
    const labelList& nPatchFaces(const label leveli) const
    {
        return nPatchFaces_[leveli];
    }

    //- Restrict (sum) a cell field with explicit addressing.
    template<class Type>
    void restrictField
    (
        Field<Type>& coarseField,
        const Field<Type>& fineField,
        const labelList& fineToCoarse
    ) const;

    //- Restrict (sum) a cell field from a fine level to the next coarse level.
    template<class Type>
    void restrictField
    (
        Field<Type>& coarseField,
        const Field<Type>& fineField,
        const label fineLevelIndex,
        const bool procAgglom = false
    ) const;

    //- Restrict (sum) a face field from a fine level to the next coarse level.
    template<class Type>
    void restrictFaceField
    (
        Field<Type>& coarseField,
        const Field<Type>& fineField,
        const label fineLevelIndex
    ) const;

    //- Prolong (inject) a cell field from a coarse level to the next fine level.
    template<class Type>
    void prolongField
    (
        Field<Type>& fineField,
        const Field<Type>& coarseField,
        const label fineLevelIndex,
        const bool procAgglom = false
    ) const;

    //- Prolong (inject) a cell field and return the coarse field used.
    template<class Type>
    const Field<Type>& prolongField
    (
        Field<Type>& fineField,
        Field<Type>& work,
        const Field<Type>& coarseField,
        const label fineLevelIndex
    ) const;


    // Processor agglomeration access

        //- Return true if processor agglomeration is enabled.
        bool doProcessorAgglomerate() const
        {
            return doProcessorAgglomerate_;
        }

        //- Return the number of coarse cells at a level.
        label nCells(const label leveli) const
        {
            return nCells_[leveli];
        }

        //- Return the number of coarse faces at a level.
        label nFaces(const label leveli) const
        {
            return nFaces_[leveli];
        }

        //- Processor restriction map at a level.
        const labelList& procAgglomMap(const label leveli) const
        {
            return procAgglomMap_[leveli];
        }

        //- Processors agglomerated onto the master at a level.
        const labelList& agglomProcIDs(const label leveli) const
        {
            return agglomProcIDs_[leveli];
        }

        //- Return true if an agglomerated mesh exists at a level.
        bool hasProcMesh(const label leveli) const
        {
            return procCommunicator_[leveli] != -1;
        }

        //- Communicator of the agglomerated mesh at a level.
        label procCommunicator(const label leveli) const
        {
            return procCommunicator_[leveli];
        }

        //- Communicator spanning the agglomerated processors at a level.
        label agglomCommunicator(const label leveli) const
        {
            return procAgglomCommunicator_[leveli].comm();
        }

        //- Cell offsets per agglomerated processor at a level.
        const labelList& cellOffsets(const label leveli) const
        {
            return procCellOffsets_[leveli];
        }

        //- Internal-face map per agglomerated processor at a level.
        const labelListList& faceMap(const label leveli) const
        {
            return procFaceMap_[leveli];
        }

        //- Boundary map per agglomerated processor at a level.
        const labelListList& boundaryMap(const label leveli) const
        {
            return procBoundaryMap_[leveli];
        }

        //- Boundary-face map per agglomerated processor at a level.
        const labelListListList& boundaryFaceMap(const label leveli) const
        {
            return procBoundaryFaceMap_[leveli];
        }


    virtual const word& type() const = 0;
};

} // End namespace Foam

#ifdef NoRepository
    #include "LduGAMGAgglomerationTemplates.C"
#endif

#endif

// ************************************************************************* //