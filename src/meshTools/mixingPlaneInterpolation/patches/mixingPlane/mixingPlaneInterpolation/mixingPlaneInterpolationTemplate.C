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

Description
    Mixing plane class dealing with transfer of data between two
    primitivePatches

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

//#include "mixingPlaneInterpolationTemplate.H"
#include "demandDrivenData.H"
#include "IOmanip.H"
#include "transform.H"
#include "vtkSurfaceWriter.H"
#include "debug.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const unsigned int
MixingPlaneInterpolation<MasterPatch, SlavePatch>::mixingPlaneVtkFeatureEdgesResolution_
(
    debug::optimisationSwitch("mixingPlaneVtkFeatureEdgesResolution", 128)
);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
Foam::direction
MixingPlaneInterpolation<MasterPatch, SlavePatch>::sweepAxisSwitch() const
{
    // SweepAxis switch
    switch (sweepAxisType_)
    {
        case SWEEP_X:
        case SWEEP_R:
        {
            return vector::X;
        }
        break;

        case SWEEP_Y:
        case SWEEP_THETA:
        {
            return vector::Y;
        }
        break;

        case SWEEP_Z:
        {
            return vector::Z;
        }
        break;

        default:
        {
            FatalErrorIn
            (
                "direction MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::sweepAxisSwitch() const"
            )   << "Bad sweepAxis type: "
                << MixingPlaneInterpolationName::sweepAxisNames_
                       [sweepAxisType_]
                << "Available types: "
                << MixingPlaneInterpolationName::sweepAxisNames_
                << abort(FatalError);

            // Dummy return
            return vector::X;
        }
    }
}


template<class MasterPatch, class SlavePatch>
Foam::direction
MixingPlaneInterpolation<MasterPatch, SlavePatch>::stackAxisSwitch() const
{
    // stackAxis switch
    switch (stackAxisType_)
    {
        case STACK_X:
        case STACK_R:
        {
            return vector::X;
        }
        break;

        case STACK_Y:
        case STACK_THETA:
        {
            return vector::Y;
        }
        break;

        case STACK_Z:
        {
            return vector::Z;
        }
        break;

        default:
        {
            FatalErrorIn
            (
                "direction MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::stackAxisSwitch() const"
            )   << "Bad stackAxis type: "
                << MixingPlaneInterpolationName::stackAxisNames_
                       [stackAxisType_]
                << "Available types: "
                << MixingPlaneInterpolationName::stackAxisNames_
                << abort(FatalError);

            // Dummy return
            return vector::X;
        }
    }
}


template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::clearOut()
{
    clearTransfomedPatches();
    clearMixingPlanePatch();

    clearAddressing();
    clearTransforms();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
MixingPlaneInterpolation
(
    const MasterPatch& masterPatch,
    const SlavePatch& slavePatch,
    const fe::coordinateSystem& cs,
    const ribbonDiscretization& discretizationType,
    const MixingPlaneInterpolationName::sweepAxis& sweepAxisType,
    const MixingPlaneInterpolationName::stackAxis& stackAxisType,
    const MixingPlaneInterpolationName::wfEvalMethod& wfEvalMethod,
    const scalar& tolPointDistanceFromRzero
)
:
    masterPatch_(masterPatch),
    slavePatch_(slavePatch),
    cs_(cs),
    ribbonDiscretization_(discretizationType),
    sweepAxisType_(sweepAxisType),
    stackAxisType_(stackAxisType),
    wfEvalMethod_(wfEvalMethod),
    tolPointDistanceFromRzero_(tolPointDistanceFromRzero),
    interpolationProfile_(),

    forwardT_(),
    reverseT_(),
    forwardSep_(),

    transformedMasterPatchPtr_(nullptr),
    transformedMasterPatchPointsPtr_(nullptr),
    transformedShadowPatchPtr_(nullptr),
    transformedShadowPatchPointsPtr_(nullptr),
    mixingPlanePatchPtr_(nullptr),
    mixingPlanePatchPointsPtr_(nullptr),

    masterPatchToProfileTPtr_(nullptr),
    masterProfileToPatchTPtr_(nullptr),
    slavePatchToProfileTPtr_(nullptr),
    slaveProfileToPatchTPtr_(nullptr),

    masterPatchToProfileAddrPtr_(nullptr),
    masterProfileToPatchAddrPtr_(nullptr),
    masterPatchToProfileWeightsPtr_(nullptr),
    masterProfileToPatchWeightsPtr_(nullptr),

    slavePatchToProfileAddrPtr_(nullptr),
    slaveProfileToPatchAddrPtr_(nullptr),
    slavePatchToProfileWeightsPtr_(nullptr),
    slaveProfileToPatchWeightsPtr_(nullptr)

{
    // Check size of transform.  They should be equal to slave patch size
    // if the transform is not constant
    if (forwardT_.size() > 1 || reverseT_.size() > 1)
    {
        if
        (
            forwardT_.size() != slavePatch_.size()
         || reverseT_.size() != masterPatch_.size()
        )
        {
            FatalErrorIn
            (
                "MixingPlaneInterpolation<MasterPatch, "
                "SlavePatch>::MixingPlaneInterpolation"
            )   << "Incorrectly defined transform: forwardT: "
                << forwardT_.size() << " patch: " << slavePatch_.size()
                << " reverseT: " << reverseT_.size()
                << " patch: " << masterPatch_.size()
                << abort(FatalError);
        }
    }

    masterPatchToProfileAddr();
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::~MixingPlaneInterpolation()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::fe::coordinateSystem&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::cs() const
{
    return cs_;
}

template<class MasterPatch, class SlavePatch>
const Foam::standAlonePatch&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::mixingPlanePatch() const
{
    if (!mixingPlanePatchPtr_)
    {
        calcMixingPlanePatch();
    }

    return *mixingPlanePatchPtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::standAlonePatch&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
transformedMasterPatch() const
{
    if (!transformedMasterPatchPtr_)
    {
        calcTransformedPatches();
    }

    return *transformedMasterPatchPtr_;
}


template<class MasterPatch, class SlavePatch>
const Foam::standAlonePatch&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
transformedShadowPatch() const
{
    if (!transformedShadowPatchPtr_)
    {
        calcTransformedPatches();
    }

    return *transformedShadowPatchPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterPatchToProfileAddr() const
{
    if (!masterPatchToProfileAddrPtr_)
    {
        calcAddressing();
    }

    return *masterPatchToProfileAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterProfileToPatchAddr() const
{
    if (!masterProfileToPatchAddrPtr_)
    {
        calcAddressing();
    }

    return *masterProfileToPatchAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterPatchToProfileWeights() const
{
    if (!masterPatchToProfileWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterPatchToProfileWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterProfileToPatchWeights() const
{
    if (!masterProfileToPatchWeightsPtr_)
    {
        calcAddressing();
    }

    return *masterProfileToPatchWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slavePatchToProfileAddr() const
{
    if (!slavePatchToProfileAddrPtr_)
    {
        calcAddressing();
    }

    return *slavePatchToProfileAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const labelListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slaveProfileToPatchAddr() const
{
    if (!slaveProfileToPatchAddrPtr_)
    {
        calcAddressing();
    }

    return *slaveProfileToPatchAddrPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slavePatchToProfileWeights() const
{
    if (!slavePatchToProfileWeightsPtr_)
    {
        calcAddressing();
    }

    return *slavePatchToProfileWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const scalarListList&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
slaveProfileToPatchWeights() const
{
    if (!slaveProfileToPatchWeightsPtr_)
    {
        calcAddressing();
    }

    return *slaveProfileToPatchWeightsPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterPatchToProfileT() const
{
    if (!masterPatchToProfileTPtr_)
    {
        calcTransforms();
    }

    return *masterPatchToProfileTPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::
masterProfileToPatchT() const
{
    if (!masterProfileToPatchTPtr_)
    {
        calcTransforms();
    }

    return *masterProfileToPatchTPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slavePatchToProfileT() const
{
    if (!slavePatchToProfileTPtr_)
    {
        calcTransforms();
    }

    return *slavePatchToProfileTPtr_;
}


template<class MasterPatch, class SlavePatch>
const tensorField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::slaveProfileToPatchT() const
{
    if (!slaveProfileToPatchTPtr_)
    {
        calcTransforms();
    }

    return *slaveProfileToPatchTPtr_;
}

template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::writePatchAsVTK
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath
) const
{
    surfaceWriters::vtkWriter writer(points, faces, outputPath);
    writer.write();
}

template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::writePointFieldAsVTK
(
    const string fName,
    const string comment,
    const pointField& points
) const
{
    // Write points into VTK file
    OFstream mps(fName + ".vtk");

    mps << "# vtk DataFile Version 4.2" << nl
        << comment << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << points.size() << " double" << nl;

    forAll (points, i)
    {
        mps << points[i].x() << ' ';
        mps << points[i].y() << ' ';
        mps << points[i].z() << ' ';

        if (i > 0 && (i % 4) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;

    return;
}

template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::writePointFieldAsCSV
(
    const string fName,
    const string comment,
    const pointField& points
) const
{
    // Write points into CSV file
    OFstream mps(fName + ".csv");

    mps << "# CSV DataFile" << nl
        << "# " << comment << nl
        << "# Number of points: " << points.size() << nl;

    forAll (points, i)
    {
        mps << points[i].x() << ',';
        mps << points[i].y() << ',';
        mps << points[i].z() << nl;
    }

    return;
}

template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::writeProfileAsVTK
(
    const string fName,
    const string comment,
    const pointField& profilePoints
) const
{
    // Write points into VTK file
    OFstream mps(fName + ".vtk");

    mps << "# vtk DataFile Version 4.2" << nl
        << comment << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << profilePoints.size() << " double" << nl;

    forAll (profilePoints, i)
    {
        mps << profilePoints[i].x() << ' ';
        mps << profilePoints[i].y() << ' ';
        mps << profilePoints[i].z() << ' ';

        if (i > 0 && (i % 4) == 0)
        {
            mps << nl;
        }
    }
    mps << nl;

    mps << "LINES 1 " << profilePoints.size() + 1 << nl;
    mps << profilePoints.size() << " ";

    forAll (profilePoints, i)
    {
        mps << i << "  ";
    }
    mps << nl;

    return;
}

template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::writeRibbonsFeatureEdgesAsVTK
(
    const string fName,
    const string comment,
    const pointField& profilePoints
) const
{
    // Get the curves resolution from tolerance switch
    const label featureEdgesRes = label(mixingPlaneVtkFeatureEdgesResolution_);

    const direction sweepDir = sweepAxisSwitch();
    const int nbrProfilePts = profilePoints.size();
    scalar sweepMin(0);
    scalar sweepMax(0);
    scalar sweepDelta(0);

    // Find the maximal possible extent of the ribbon patch in the
    // sweepDir direction
    if (cs().spanLimited()[sweepDir][0] && cs().spanLimited()[sweepDir][1])
    {
        // Span limited along the sweepDir: use the coordinates system
        // spanBounds values
        sweepMin = cs().spanBounds().min()[sweepDir];
        sweepMax = cs().spanBounds().max()[sweepDir];
    }
    else
    {
        // Not span limited: use the mixingPlane patches bounding box along the
        // sweepDir direction in local coords
        boundBox masterBB(cs().localPosition(masterPatch_.localPoints()), false);
        boundBox shadowBB(cs().localPosition(slavePatch_.localPoints()), false);
        boundBox globSpanBB
        (
            point(Foam::min(masterBB.min(), shadowBB.min())),
            point(Foam::max(masterBB.max(), shadowBB.max()))
        );

        sweepMin = globSpanBB.min()[sweepDir];
        sweepMax = globSpanBB.max()[sweepDir];
    }

    sweepDelta = (sweepMax - sweepMin) / featureEdgesRes;

    // Write curves into VTK file, one curve per profile point
    OFstream mps(fName + ".vtk");

    mps << "# vtk DataFile Version 4.2" << nl
        << comment << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << nbrProfilePts * (featureEdgesRes + 1)
        << " double" << nl;

    // Generate the feature edges in global coordinates, located at each profile
    // point
    for (auto localPoint : profilePoints)
    {
        scalar sweepAngle = sweepMin;

        for (int sweepI=0; sweepI<=featureEdgesRes; sweepI++)
        {
            localPoint[sweepDir] = sweepAngle;
            point globalPoint = cs().globalPosition(localPoint);

            mps << globalPoint.x() << ' ';
            mps << globalPoint.y() << ' ';
            mps << globalPoint.z() << ' ';

            if (sweepI > 0 && (sweepI % 5) == 0)
            {
                mps << nl;
            }

            sweepAngle += sweepDelta;
        }
    }
    mps << nl;

    mps << "LINES "  << nbrProfilePts << " "
        << nbrProfilePts * (featureEdgesRes + 2)
        << nl;

    for (int i=0; i< nbrProfilePts; i++)
    {
        mps << featureEdgesRes + 1 << " ";

        for (int sweepI=0; sweepI<=featureEdgesRes; sweepI++)
        {
            mps << (featureEdgesRes+1)*i + sweepI << " ";
        }
        mps << nl;
    }
    mps << nl;

    return;
}

template<class MasterPatch, class SlavePatch>
void MixingPlaneInterpolation<MasterPatch, SlavePatch>::writeMixingPlaneRibbonGenerationData
(
    const fileName& caseDirName,
    const word& masterPatchName,
    const word& neighbourPatchName
) const
{
    // Create dump directories for VTK and CSV files
    fileName vtkDirName(caseDirName/"VTK");
    fileName csvDirName(caseDirName/"CSV");

    mkDir(vtkDirName);
    mkDir(csvDirName);

    const string interfaceIdForFileName
    (
        "mixingPlane_" +
        masterPatchName +
        "_" +
        neighbourPatchName
    );
    const word interfaceIdForConsoleMsg
    (
        " Master: " +
        masterPatchName +
        " Neighbour: " +
        neighbourPatchName
    );

    Info<< "Writing mixing plane patches in local coordinates as VTK."
        << interfaceIdForConsoleMsg << endl;

    // MixingPlane patch
    writePatchAsVTK
    (
        mixingPlanePatch().localPoints(),
        mixingPlanePatch().localFaces(),
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_ribbonPatch"
        )
    );

    // Master patch
    writePatchAsVTK
    (
        transformedMasterPatch().localPoints(),
        transformedMasterPatch().localFaces(),
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_masterPatch"
        )
    );

    // Neighbour patch
    writePatchAsVTK
    (
        transformedShadowPatch().localPoints(),
        transformedShadowPatch().localFaces(),
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_neighbourPatch"
        )
    );

    // Write the ribbon patch profile information in both CSV and VTK format
    // The CSV files can be used for developping user-defined ribbons patch
    // profiles
    pointField masterGlobalProfile
    (
        cs().localPosition(masterPatch_.localPoints())
    );

    pointField shadowGlobalProfile
    (
        cs().localPosition(slavePatch_.localPoints())
    );

    const direction sweepDir = sweepAxisSwitch();

    // Collapse all points in the sweep direction
    masterGlobalProfile.replace(sweepDir, 0);
    shadowGlobalProfile.replace(sweepDir, 0);

    Info << "Writing master and neighbour raw profile points as VTK in global coordinates."
        << interfaceIdForConsoleMsg << endl;

    writePointFieldAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_masterPatch_rawProfilePoints"
        ),
        "mixingPlane raw profile points in global coordinates: from master patch",
        cs().globalPosition(masterGlobalProfile)
    );
    writePointFieldAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_neighbourPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in global coordinates: from neighbour patch",
        cs().globalPosition(shadowGlobalProfile)
    );

    Info << "Writing master and neighbour raw profile points as CSV in global coordinates."
        << interfaceIdForConsoleMsg << endl;

    writePointFieldAsCSV
    (
        csvDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_masterPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in global coordinates: from master patch",
        cs().globalPosition(masterGlobalProfile)
    );
    writePointFieldAsCSV
    (
        csvDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_neighbourPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in global coordinates: from neighbour patch",
        cs().globalPosition(shadowGlobalProfile)
    );

    Info << "Writing master and neighbour raw profile points as VTK in local coordinates."
        << interfaceIdForConsoleMsg << endl;

    writePointFieldAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_masterPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in local coordinates: from master patch",
        masterGlobalProfile
    );
    writePointFieldAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_neighbourPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in local coordinates: from neighbour patch",
        shadowGlobalProfile
    );

    Info << "Writing master and neighbour raw profile points as CSV in local coordinates."
        << interfaceIdForConsoleMsg << endl;

    writePointFieldAsCSV
    (
        csvDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_masterPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in local coordinates: from master patch",
        masterGlobalProfile
    );
    writePointFieldAsCSV
    (
        csvDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_neighbourPatch_rawProfilePoints"
        ),
        "mixingPlane profile raw points in local coordinates: from neighbour patch",
        shadowGlobalProfile
    );

    Info << "Writing mixingPlane profile as VTK in global coordinates."
        << interfaceIdForConsoleMsg << endl;

    writeProfileAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_profile"
        ),
        "mixingPlane profile as polyline in global coordinates",
        cs().globalPosition(interpolationProfile_)
    );

    Info << "Writing mixingPlane profile as CSV in global coordinates."
        << interfaceIdForConsoleMsg << endl;

    writePointFieldAsCSV
    (
        csvDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_profile"
        ),
        "mixingPlane profile as polyline in global coordinates",
        cs().globalPosition(interpolationProfile_)
    );

    Info << "Writing mixingPlane profile as VTK in local coordinates."
        << interfaceIdForConsoleMsg << endl;

    writeProfileAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_profile"
        ),
        "mixingPlane profile as polyline in local coordinates",
        interpolationProfile_
    );

    Info << "Writing mixingPlane profile as CSV in local coordinates."
        << interfaceIdForConsoleMsg << endl;

    writePointFieldAsCSV
    (
        csvDirName/fileName
        (
            interfaceIdForFileName + "_localCoords_profile"
        ),
        "mixingPlane profile as polyline in local coordinates",
        interpolationProfile_
    );

    // Write ribbons patch feature edges as a VTK file
    Info << "Writing mixingPlane ribbons patch feature edges as VTK."
        << interfaceIdForConsoleMsg << endl;

    writeRibbonsFeatureEdgesAsVTK
    (
        vtkDirName/fileName
        (
            interfaceIdForFileName + "_globalCoords_ribbonPatchFeatureEdges"
        ),
        "mixingPlane ribbon patch feature edges",
        interpolationProfile_
    );

    return;
}

template<class MasterPatch, class SlavePatch>
bool MixingPlaneInterpolation<MasterPatch, SlavePatch>::movePoints()
{
    clearOut();

    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
