/*---------------------------------------------------------------------------*\

License
    This file is part of the TurboWG project:
    https://sourceforge.net/projects/turbowg

    The files from the TurboWG project are free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    This source code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this source code.  If not, see <http://www.gnu.org/licenses/>.

    Special notice:
    This file is part of a port of the mixingPlane interface source code from
    the foam-extend project (https://sourceforge.net/projects/foam-extend/) to
    OpenFOAM (https://www.openfoam.com/). Most of the original source code has
    been preserved, including the comments and authorship/IP notices.

Contributors/Copyright:
    Martin Beaudoin, (2023)

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     5.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    mixingPlanePolyPatch

Description
    Mixing plane patch

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlanePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlanePolyPatch_H
#define mixingPlanePolyPatch_H

#include "coupledPolyPatch.H"
#include "mixingPlaneInterpolation.H"
#include "transform.H"
#include "word.H"
#include "faceZone.H"
#include "fe_coordinateSystem.H"
#include "coupleGroupIdentifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mixingPlanePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingPlanePolyPatch
:
    public coupledPolyPatch
{
    // Private data

         //- Name of neighbour patch
        mutable word neighbPatchName_;

        //- Shadow patch name
        const word shadowName_;

        //- Interpolation zone name
        const word zoneName_;

        //- Index of other half
        mutable label neighbPatchID_;

        //- Optional patchGroup to find neighbPatch
        const coupleGroupIdentifier coupleGroup_;

        //- Coordinate system for the mixing plane
        autoPtr<fe::coordinateSystem> csPtr_;

        //- Type of mixing plane discretization algorithm
        ribbonDiscretization ribbonDiscretization_;

        //- Orientation of the mixing plane ribbon patch
        //  sweepAxis = axis in which ribbon length is measured
        //  stackAxis = axis in which ribbon width/height is measured
        mixingPlaneInterpolation::sweepAxis sweepAxisType_;
        mixingPlaneInterpolation::stackAxis stackAxisType_;

        //- Method for evaluating the mixingPlane weighting factors
        //  GGI (default), AMI
        mixingPlaneInterpolation::wfEvalMethod wfEvalMethod_;

        //- User-defined interpolation profile file name
        fileName userProfileFile_;

        //- When saving mixingPlane patchFields to file, force the update of the
        //  mixingPlane patchField before writing them to disk.
        //  This boolean flag can be specified for either the master or shadow
        //  mixingPlane patches.
        //  If left unspecified, or set to false, the default behaviour is to
        //  write the patchField values as they were left after the solver's latest
        //  iteration; those patchField values might no longer be axisymmetric.
        //  If set to true, the "areaAveraged" patchFields will be re-evaluated
        //  and will appear as axisymmetric, while the "fluxAveraging" patchFields
        //  will be re-evaluated according to the latest flux values/directions.
        //  This option is useful mainly for visualization and illustration purposes.
        bool writeAveragedPatchField_;

        // Demand-driven data

            //- Shadow patch index.  Delayed evaluation for construction
            mutable label shadowIndex_;

            //- Interpolation zone index.  Delayed evaluation for construction
            mutable label zoneIndex_;

            //- Patch-to-patch interpolation
            mutable mixingPlaneZoneInterpolation* patchToPatchPtr_;

            //- Patch-to-zone addressing
            mutable labelList* zoneAddressingPtr_;

            //- Reconstructed patch neighbour cell centres
            mutable vectorField* reconFaceCellCentresPtr_;


        // Parallel communication optimisation, stored on master processor

            //- Is the patch localised on a single processor
            // (single processor in a parallel run)?
            //  Used for parallel optimisation
            mutable bool* localParallelPtr_;

            //- List of zone faces indices received from each processor
            mutable labelListList* receiveAddrPtr_;

    // Private member functions

        //- Calculate patch-to-zone addressing
        virtual void calcZoneAddressing() const;

        //- Calculate interpolation
        virtual void calcPatchToPatch() const;


        // Geometry

            //- Calculate reconstructed cell centres
            void calcReconFaceCellCentres() const;

            //- Force calculation of transformation tensors
            virtual void calcTransforms() const;


        // Parallel communication optimisation, stored on master processor

            //- Calculate local parallel switch
            void calcLocalParallel() const;

            //- Calculate receive addressing
            void calcReceive() const;


            //- Return receive addressing
            const labelListList& receiveAddr() const;


        // Memory management

            //- Clear geometry
             void clearLocalGeom() const;

            //- Clear out
            void clearOut() const;


protected:

    // Protected Member functions

        //- Is the mixing plane active? (zone and shadow present)
        bool active() const;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry(PstreamBuffers&);

        //- Calculate the patch geometry
        virtual void calcGeometry
        (
            const primitivePatch& referPatch,
            const pointField& thisCtrs,
            const vectorField& thisAreas,
            const pointField& thisCc,
            const pointField& nbrCtrs,
            const vectorField& nbrAreas,
            const pointField& nbrCc
        );

        //- Initialise the patches for moving points
        virtual void initMovePoints(PstreamBuffers&, const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(PstreamBuffers&, const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh(PstreamBuffers&);

        //- Update of the patch topology
        virtual void updateMesh(PstreamBuffers&);


public:

    //- Runtime type information
    TypeName("mixingPlane");

    // Static data

        //- Tolerance distance for evaluating facet points proximity to r=0
        static scalar tolPointDistanceFromRzero_;

        //- Info switch: write ribbons profile information to file (CSV and VTK)
        static bool writeRibbonsProfileToFile_;


    // Constructors

        //- Construct from components
        mixingPlanePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType,
            const transformType transform = UNKNOWN
        );

        //- Construct from dictionary
        mixingPlanePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm,
            const word& patchType
        );

        //- Construct as copy, resetting the boundary mesh
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm
        );

        //- Copy construct, resetting nbrPatchID and faceCells
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch& pp,
            label nbrPatchID,
            const labelList& faceCells
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart,
            const word& neighbPatchName
        );

        //- Construct given the original patch and a map
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new mixingPlanePolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new mixingPlanePolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart,
                    neighbPatchName_
                )
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new mixingPlanePolyPatch
                (
                    *this,
                    bm,
                    index,
                    mapAddressing,
                    newStart
                )
            );
        }

    //- Destructor
    virtual ~mixingPlanePolyPatch();


    // Member functions

            //- Is patch 'coupled'. Note that with mixingPlane, like the AMI interface,
            //- the geometry is not coupled but the fields are!
            virtual bool coupled() const
            {
                return false;
            }

        // Basic info

            //- Return shadow patch name
            const word& shadowName() const
            {
                return shadowName_;
            }

            //- Return name of interpolation face zone
            const word& zoneName() const
            {
                return zoneName_;
            }

            //- Return shadow patch index
            label shadowIndex() const;

            //- Return zone patch index
            label zoneIndex() const;

            //- Return shadow patch
            const mixingPlanePolyPatch& shadow() const;

            //- Return interpolation face zone
            const faceZone& zone() const;

            //- Transform a patch-based position from other side to this side
            virtual void transformPosition(pointField& l) const;

            //- Transform a patch-based position from other side to this side
            virtual void transformPosition(point&, const label facei) const;

            //- Return the writeAveragedPatchField_ flag value
            bool writeAveragedPatchField() const
            {
	        return writeAveragedPatchField_;
            }


        // Interpolation data

            //- Is this the master side?
            virtual bool master() const
            {
                return index() < shadowIndex();
            }

            //- Neighbour patch name
            const word& neighbPatchName() const;

            //- Neighbour patchID
            virtual label neighbPatchID() const;

            virtual bool owner() const
            {
                return index() < neighbPatchID();
            }

            virtual bool neighbour() const
            {
                return !owner();
            }

            //- Return coordinate system
            const fe::coordinateSystem& cs() const;

            //- Return number of profile bands
            label nProfileBands() const
            {
                return patchToPatch().nProfileBands();
            }

            //- Return patch-to-zone addressing
            const labelList& zoneAddressing() const;

            //- Is the patch localised on a single processor
            bool localParallel() const;

            //- Return reference to patch-to-patch interpolation
            //  Used only for addressing
            const mixingPlaneZoneInterpolation& patchToPatch() const;


        // Communications support

            //- Return communicator used for communication
            virtual label comm() const;

            //- Return message tag to use for communication
            virtual int tag() const;


        // Interpolation functions

            //- Interpolate face field to profile: given field on a
            //- master/slave side, create a profile field
            template<class Type>
            tmp<Field<Type> > toProfile
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > toProfile
            (
                const tmp<Field<Type> >& tpf
            ) const;

            //- Interpolate profile to face field: given field on a profile,
            //  create a patch field
            template<class Type>
            tmp<Field<Type> > fromProfile
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > fromProfile
            (
                const tmp<Field<Type> >& tpf
            ) const;

            //- Interpolate face field: given field on a the shadow side,
            //  create an interpolated field on this side
            template<class Type>
            tmp<Field<Type> > interpolate(const Field<Type>& pf) const;

            template<class Type>
            tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const;

            //- Interpolate face field from average: given field on a the
            //  master/slave side, create interpolated field on the same side
            //  using averaged values
            template<class Type>
            tmp<Field<Type> > circumferentialAverage
            (
                const Field<Type>& pf
            ) const;

            template<class Type>
            tmp<Field<Type> > circumferentialAverage
            (
                const tmp<Field<Type> >& tpf
            ) const;


        // Geometric data

            //- Return reconstructed cell centres
            const vectorField& reconFaceCellCentres() const;


        // Patch ordering

            //- Initialize ordering for primitivePatch. Does not
            //  refer to *this (except for name() and type() etc.)
            virtual void initOrder(PstreamBuffers&, const primitivePatch&) const;

            //- Return new ordering for primitivePatch.
            //  Ordering is -faceMap: for every face
            //  index of the new face -rotation: for every new face the
            //   clockwise shift of the original face.
            //   Return false if nothing changes
            //  (faceMap is identity, rotation is 0), true otherwise.
            virtual bool order
            (
                PstreamBuffers&,
                const primitivePatch&,
                labelList& faceMap,
                labelList& rotation
            ) const;

            //- Synchronise communications of ordering for primitivePatch
            //  Used in cases when no topological change happens locally,
            //  but is happening on other processors
            virtual void syncOrder() const;

        //- Check functions

        // Used to validate and report about specific geometry and topology
        // configurations for coupled regions, like presence of non-overlapped
        // GGI facets, misaligned rotation axis for cyclicGGI patches, or
        // absence of a necessary discretization point at r=0 for some
        // configurations of the master/shadow patches for the mixingPlane
        // interface.
        // Will be called by checkMesh with report flag enabled.

        // Check for geometry problems
        virtual bool checkGeometry(const bool report=false) const;

        // Check for topology problems
        virtual bool checkTopology(const bool report=false) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingPlanePolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
