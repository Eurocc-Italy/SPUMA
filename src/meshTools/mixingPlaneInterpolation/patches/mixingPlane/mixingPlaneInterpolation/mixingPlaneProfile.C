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
    MixingPlane interpolation functions

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
tmp<pointField>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::computeProfileFromHistograms
(
    const profileHistogram& masterHisto,
    const profileHistogram& slaveHisto,
    const scalar halfSizeBin
) const
{
    // Find min, max bounds
    scalar histoMinValue =
        Foam::min
        (
            masterHisto.begin()->first,
            slaveHisto.begin()->first
        );

    scalar histoMaxValue =
        Foam::max
        (
            (--masterHisto.end())->first,
            (--slaveHisto.end())->first
        );

    // Next, we compare both histograms, leap-frogging from one histogram
    // to the other, everytime jumping to the next largest value from
    // a given position
    std::list<point> leapFrogProfile(0);

    scalar curRvalue = histoMinValue;

    if (masterHisto.find(curRvalue) != masterHisto.end())
    {
        leapFrogProfile.push_back
        (
            *(masterHisto.find(curRvalue)->second.begin())
        );
    }
    else
    {
        leapFrogProfile.push_back
         (
             *(slaveHisto.find(curRvalue)->second.begin())
         );
    }

    do
    {
        profileHistogram::const_iterator nextMaster =
            masterHisto.lower_bound(curRvalue + halfSizeBin - SMALL);

        profileHistogram::const_iterator nextSlave  =
            slaveHisto.lower_bound(curRvalue + halfSizeBin - SMALL);

        if
        (
            nextMaster == masterHisto.end()
         || nextSlave  == slaveHisto.end()
        )
        {
            // We are done
            if (curRvalue != histoMaxValue)
            {
                if (masterHisto.find(histoMaxValue) != masterHisto.end())
                {
                    leapFrogProfile.push_back
                    (
                        *(masterHisto.find(histoMaxValue)->second.begin())
                    );
                }
                else
                {
                    leapFrogProfile.push_back
                    (
                        *(slaveHisto.find(histoMaxValue)->second.begin())
                    );
                }
            }
            break;
        }

        // Leap frog to the next largest delta
        // Please note that we are not taking into account the number of
        // values at the specific map index (the attribute called 'second'
        // of the map container)
        // So by leap-frogging this way, we might ran into the situation
        // where the next jump will bring us to a radius value shared by only
        // one patch face point, and we might oversee a slightly smaller
        // radius shared by 100s of points.
        // There is a definite advantage of sticking to existing radius values
        // for the profile because this will tend to keep the number of
        // GGI face neighbours low.
        // I am not sure if there is a nice solution to this problem.... MB
        //
        // Also note also that the profile extremum values are stored as the
        // first entry in their respective bin, which are the first and last
        // bins of the histograms. MB 03/2021
        //
        curRvalue = max(nextMaster->first, nextSlave->first);

        if (masterHisto.find(curRvalue) != masterHisto.end())
        {
            leapFrogProfile.push_back
            (
                *(masterHisto.find(curRvalue)->second.begin())
            );
        }
        else
        {
            leapFrogProfile.push_back
            (
                *(slaveHisto.find(curRvalue)->second.begin())
            );
        }
    } while (curRvalue != histoMaxValue);

    // Re-package the data into pointField
    tmp<pointField> tprofile(new pointField(leapFrogProfile.size()));

    pointField& profile = tprofile.ref();

    label pI = 0;

    forAllIter (std::list<point>, leapFrogProfile, lI)
    {
        profile[pI++] = *lI;
    }

    return tprofile;
}


template<class MasterPatch, class SlavePatch>
void
MixingPlaneInterpolation<MasterPatch, SlavePatch>::createPointCoordsHistogram
(
    const pointField& profileCoords, // raw 3D point coordinates
    const direction dir,             // Sorting dimension 0: x, 1: y, 2: z
    const scalar halfSizeBin,        // half size of min width for histogram bins
    profileHistogram& histo          // histogram container to initialize
) const
{
    // Initial cleanup
    histo.clear();

    // Find the profile extrema, making sure the profile will cover the full
    // span in the dir direction

    point minPoint = profileCoords[0];
    point maxPoint = profileCoords[0];

    forAll (profileCoords, pI)
    {
        scalar keyValue = profileCoords[pI].component(dir);

        if (keyValue > maxPoint.component(dir))
        {
            maxPoint = profileCoords[pI];
        }
        else if (keyValue < minPoint.component(dir))
        {
            minPoint = profileCoords[pI];
        }
    }

    // Create a bin for each of those extrema, making sure they are the first
    // of the bin's list of points
    std::list<point> minList = { minPoint };
    std::list<point> maxList = { maxPoint };

    histo.insert
        (
            std::pair<scalar, std::list<point> >
            (
                minPoint.component(dir),
                minList
            )
        );

    histo.insert
        (
            std::pair<scalar, std::list<point> >
            (
                maxPoint.component(dir),
                maxList
            )
        );


    // Iterate over coordinates, packed into bins not larger than 2*halfSizeBin
    forAll (profileCoords, pI)
    {
        bool foundNewBin = true;

        scalar keyValue = profileCoords[pI].component(dir);

        forAllIter (profileHistogram, histo, histoI)
        {
            if
            (
                keyValue >= histoI->first - halfSizeBin
                && keyValue < histoI->first + halfSizeBin
            )
            {
                foundNewBin = false;
                histoI->second.push_back(profileCoords[pI]);
                break;
            }
        }

        if (foundNewBin)
        {
            std::list<point> initValue;

            initValue.push_back(profileCoords[pI]);
            histo.insert
            (
                std::pair<scalar, std::list<point> >(keyValue, initValue)
            );
        }
    }

    if (debug > 1)
    {
        // Dump histogram info to console
        forAllConstIter (profileHistogram, histo, zHi)
        {
            Info<< "histo (z, n): (" << zHi->first << " "
                << static_cast<int>(zHi->second.size()) << ")" << endl;
        }
    }

    return;
}

template<class MasterPatch, class SlavePatch>
tmp<pointField>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::computeLinearProfileAlongStackAxis
(
    const boundBox& globBoundingBox,
    const direction stackingDir,
    const direction sweepDir,
    const label nbrRibbons
) const
{
    // Allocate profile
    tmp<pointField> tprofile(new pointField(nbrRibbons + 1));

    pointField& profile = tprofile.ref();

    point delta(0.0, 0.0, 0.0);

    // Compute the size of each ribbon along the stacking direction
    delta[stackingDir] =
        (
            globBoundingBox.max()[stackingDir] -
            globBoundingBox.min()[stackingDir]
        ) / nbrRibbons;

    // Simply add a sequence of delta offsets from the starting profile point.
    // The generated profile will follow purely horizontal and purely vertical
    // flat master/slave patches. For strongly curved or oblique geometry, you
    // might have to adjust the runtime GGIFeatureCosTol tolerance factor in
    // order to account for larger variations between potential GGI facets
    // normals.
    // MB 03/2021

    profile[0] = globBoundingBox.min();

    for(label i=1; i<=nbrRibbons; i++)
    {
        profile[i] = profile[i-1] + delta;
    }

    // Make sure the profile coordinates are all aligned along sweepDir
    profile.replace(sweepDir, 0);

    return tprofile;
}

template<class MasterPatch, class SlavePatch>
tmp<pointField>
MixingPlaneInterpolation<MasterPatch, SlavePatch>::calcProfile() const
{
    // This is the main container for the ribbon profile. Every processor will
    // have a duplicate of this profile so they can generate an identical ribbon
    // patch on each processor.
    pointField profileBeforeClip;

    // First, transform patch points over one single global profile
    // using local coordinates

    pointField masterPatchPoints
    (
        cs().localPosition(masterPatch_.localPoints())
    );

    pointField slavePatchPoints
    (
        cs().localPosition(slavePatch_.localPoints())
    );

    // Check for the presence of a central facet with no discretization
    // coordinates at r=0.
    // If so, we propose a solution for fixing the mesh, and we abort the computation.
    // We only do this test for radius-based coordinates systems, with an exception for
    // the toroidal coordinates system.
    // MB 08/2023
    if
    (
        cs().type() != "coordinateSystem" &&
        cs().type() != "toroidal"
    )
    {
        checkMissingPatchDiscretizationAtRzero
        (
            masterPatchPoints,
            slavePatchPoints
        );
    }

    // Get stackingDir and sweepDir switches
    const direction stackingDir = stackAxisSwitch();
    const direction sweepDir = sweepAxisSwitch();

    // Gather patch points coordinates onto master processor

    List<pointField> procMasterPatchCoords(Pstream::nProcs());
    List<pointField> procSlavePatchCoords(Pstream::nProcs());
    procMasterPatchCoords[Pstream::myProcNo()] = masterPatchPoints;
    procSlavePatchCoords[Pstream::myProcNo()] = slavePatchPoints;

    // Gather data onto master processor
    Pstream::gatherList(procMasterPatchCoords);
    Pstream::gatherList(procSlavePatchCoords);


    // Next, find the smallest edge length from both patches.
    // This length will control the size of the bins for the histograms
    // we are about to build.
    scalar masterMinEdgeLength = GREAT;
    const edgeList& masterEdgeList = masterPatch_.edges();

    forAll (masterEdgeList, mEi)
    {
        masterMinEdgeLength =
            Foam::min
            (
                masterMinEdgeLength,
                masterEdgeList[mEi].mag(masterPatch_.localPoints())
            );
    }

    scalar slaveMinEdgeLength = GREAT;
    const edgeList& slaveEdgeList = slavePatch_.edges();

    forAll (slaveEdgeList, sEi)
    {
        slaveMinEdgeLength =
            Foam::min
            (
                slaveMinEdgeLength,
                slaveEdgeList[sEi].mag(slavePatch_.localPoints())
            );
    }

    // Reduce min edge values among partitions
    if(Pstream::parRun())
    {
        masterMinEdgeLength = returnReduce(masterMinEdgeLength, minOp<scalar>());
        slaveMinEdgeLength  = returnReduce(slaveMinEdgeLength, minOp<scalar>());
    }

    // There is no point classifying data to a resolution smaller than the
    // largest of the two minimum edges found.
    // This will drive the size of the smallest bin needed to construct the
    // z and r histograms
    scalar halfMinSizeBin =
        Foam::max(masterMinEdgeLength, slaveMinEdgeLength)/2.0;

    if (debug)
    {
        InfoInFunction
            << "halfMinSizeBin: " << halfMinSizeBin << endl;
    }

    // Next, evaluate the ribbon patch profile on the master processor
    if (Pstream::master())
    {
        // Get profile patch points from all processors
        pointField masterGlobalProfile = ListListOps::combine<pointField>
            (
                procMasterPatchCoords,
                accessOp<pointField>()
            );

        pointField slaveGlobalProfile = ListListOps::combine<pointField>
            (
                procSlavePatchCoords,
                accessOp<pointField>()
            );

        //Pout << "DEBUG:  masterPatchPoints.size() " << masterPatchPoints.size() << endl;
        //Pout << "DEBUG:  totalSize: " << masterGlobalProfile.size() << endl;

        // Collapse all points in the sweep direction
        masterGlobalProfile.replace(sweepDir, 0);
        slaveGlobalProfile.replace(sweepDir, 0);

        // Compute the global bounding box master/slave patch points
        boundBox masterBB(masterGlobalProfile, false);
        boundBox slaveBB (slaveGlobalProfile,  false);

        // Expand the bounding box in the sweepAxis-wise direction
        // Note: All points are collapsed to zero in sweepDir
        // It is sufficient to expand the box by 1 in this direction

        masterBB.min().replace(sweepDir, -1);
        masterBB.max().replace(sweepDir, 1);

        slaveBB.min().replace(sweepDir, -1);
        slaveBB.max().replace(sweepDir, 1);

        boundBox globSpanBB
        (
            point(Foam::min(masterBB.min(), slaveBB.min()))
            - point(SMALL, SMALL, SMALL),
            point(Foam::max(masterBB.max(), slaveBB.max()))
            + point(SMALL, SMALL, SMALL)
        );

        if (debug)
        {
            InfoInFunction
                << "masterGlobalProfile: " << masterGlobalProfile << nl
                << "slaveGlobalProfile: " << slaveGlobalProfile << endl;
        }

        // Select which type of profile we need
        switch (ribbonDiscretization_.method)
        {
            case MASTER_PATCH:   // Deprecated
            {
                WarningInFunction
                    << "Using deprecated 'discretization masterPatch' method "
                    << "instead of 'followPatchDiscretization' method"
                    << nl << endl;

                profileHistogram masterHistogram;
                createPointCoordsHistogram
                (
                    masterGlobalProfile,
                    stackingDir,
                    halfMinSizeBin,
                    masterHistogram
                );

                profileBeforeClip = computeProfileFromHistograms
                (
                    masterHistogram,
                    masterHistogram,
                    halfMinSizeBin
                );
            }
            break;

            case SLAVE_PATCH:   // Deprecated
            {
                WarningInFunction
                    << "Using deprecated 'discretization slavePatch' method "
                    << "instead of 'followPatchDiscretization' method"
                    << nl << endl;

                profileHistogram slaveHistogram;
                createPointCoordsHistogram
                (
                    slaveGlobalProfile,
                    stackingDir,
                    halfMinSizeBin,
                    slaveHistogram
                );

                profileBeforeClip = computeProfileFromHistograms
                (
                    slaveHistogram,
                    slaveHistogram,
                    halfMinSizeBin
                );
            }
            break;

            case BOTH_PATCHES:   // Deprecated
            {
                WarningInFunction
                    << "Using deprecated 'discretization bothPatches' method "
                    << "instead of 'followPatchDiscretization' method"
                    << nl << endl;

                profileHistogram masterHistogram;
                createPointCoordsHistogram
                (
                    masterGlobalProfile,
                    stackingDir,
                    halfMinSizeBin,
                    masterHistogram
                );

                profileHistogram slaveHistogram;
                createPointCoordsHistogram
                (
                    slaveGlobalProfile,
                    stackingDir,
                    halfMinSizeBin,
                    slaveHistogram
                );

                profileBeforeClip = computeProfileFromHistograms
                (
                    masterHistogram,
                    slaveHistogram,
                    halfMinSizeBin
                );
            }
            break;

            case FOLLOW_PATCH_DISCRETIZATION:
            {
                // Extract user parameters
                dictionary& dict = ribbonDiscretization_.methodParamsDict;

                // Get discretization option
                followPatchDiscretizationOption followOption =
                    MixingPlaneInterpolationName::followPatchDiscretizationOptionNames_
                    [
                        dict.lookupOrDefault<word>("follow", "Undefined")
                    ];

                // Get the optional histogram's bin size scale factor
                scalar histoBinSizeScalingFactor =
                    dict.lookupOrDefault<scalar>("histoBinSizeScaleFactor", 1.0);

                halfMinSizeBin *= histoBinSizeScalingFactor;

                switch(followOption)
                {
                    case FOLLOW_MASTER_PATCH:
                    {
                        profileHistogram masterHistogram;
                        createPointCoordsHistogram
                        (
                            masterGlobalProfile,
                            stackingDir,
                            halfMinSizeBin,
                            masterHistogram
                        );

                        profileBeforeClip = computeProfileFromHistograms
                        (
                            masterHistogram,
                            masterHistogram,
                            halfMinSizeBin
                        );
                    }
                    break;

                    case FOLLOW_NEIGHBOUR_PATCH:
                    {
                        profileHistogram slaveHistogram;
                        createPointCoordsHistogram
                        (
                            slaveGlobalProfile,
                            stackingDir,
                            halfMinSizeBin,
                            slaveHistogram
                        );

                        profileBeforeClip = computeProfileFromHistograms
                        (
                            slaveHistogram,
                            slaveHistogram,
                            halfMinSizeBin
                        );
                    }
                    break;

                    case FOLLOW_BOTH_PATCHES:
                    {
                        profileHistogram masterHistogram;
                        createPointCoordsHistogram
                        (
                            masterGlobalProfile,
                            stackingDir,
                            halfMinSizeBin,
                            masterHistogram
                        );

                        profileHistogram slaveHistogram;
                        createPointCoordsHistogram
                        (
                            slaveGlobalProfile,
                            stackingDir,
                            halfMinSizeBin,
                            slaveHistogram
                        );

                        profileBeforeClip = computeProfileFromHistograms
                        (
                            masterHistogram,
                            slaveHistogram,
                            halfMinSizeBin
                        );
                    }
                    break;

                    default:
                    {
                        FatalErrorInFunction
                            << "Bad followPatchDiscretization option: "
                            << followPatchDiscretizationOptionNames_[followOption]
                            << nl << "Available options are: " << nl
                            << followPatchDiscretizationOptionNames_
                            << abort(FatalError);
                    }
                    break;
                }

            }
            break;

            case LINEAR_ALONG_STACK_AXIS:
            {
                // Extract user parameters
                dictionary& dict = ribbonDiscretization_.methodParamsDict;

                label nRibbons = dict.lookupOrDefault<label>("nRibbons", 0);
                scalarListList grading = dict.lookupOrDefault<scalarListList>
                (
                    "grading",
                    scalarListList::null()
                );

                // Warn about grading not yet implemented
                if(grading.size() > 0)
                {
                    WarningInFunction
                        << ": grading of profile not implemented yet." << endl;
                }

                profileBeforeClip = computeLinearProfileAlongStackAxis
                (
                    globSpanBB,
                    stackingDir,
                    sweepDir,
                    nRibbons
                );

                if(debug)
                    InfoInFunction << "linear profile: " << profileBeforeClip << endl;
            }
            break;

            case USER_DEFINED:
            {
                // Extract user parameters
                dictionary& dict = ribbonDiscretization_.methodParamsDict;

                fileName userProfileFile(dict.lookup("file"));

                Info << "Reading mixingPlane ribbons profile from file: "
                     << userProfileFile << endl
                     << endl;

                // Construct CSV file reading object from dictionary
                Function1Types::twg_CSV<Foam::point> csvProfile("profile", dict);

                // Grab the user-defined profile coordinates
                vectorField csvData(csvProfile.y());

                if (csvData.empty())
                {
                    FatalErrorInFunction
                        << "Empty user-defined mixing plane profile read from file "
                        << userProfileFile
                        << abort(FatalError);
                }

                // Check for interpolation of profile
                label nbrRequestedRibbons = dict.lookupOrDefault("nRibbons", 0);

                if(nbrRequestedRibbons > 0)
                {
                    // Not implemented yet.
                    // Was caught earlier in mixingPlanePolyPatch
                }

                // Check for the global coordinate system flag. The profile
                // coordinates must be converted in the local coordinate system
                userDefinedProfilecoordinateSystemType coordsType =
                    userDefinedProfilecoordinateSystemTypeNames_
                    [
                        dict.lookupOrDefault<word>("coordsType", "global")
                    ];

                if(coordsType == USER_DEFINED_PROFILE_COORDS_GLOBAL)
                {
                    // Convert coords to local coordinates system
                    profileBeforeClip = cs().localPosition(csvData);
                }
                else
                {
                    // Keep as is
                    profileBeforeClip = csvData;
                }

                // Make sure the profile coordinates are all aligned along sweepDir
                profileBeforeClip.replace(sweepDir, 0);

                if(debug)
                {
                    InfoInFunction << "User-defined profile from file: "
                        << userProfileFile
                        << " : " << nl
                        << csvData << endl;

                    InfoInFunction << "User-defined profile in local coordinates: "
                        << profileBeforeClip << endl;
                }
            }
            break;

            case POLYNOMIAL:
            {
                notImplemented
                (
                    string(FUNCTION_NAME) + nl +
                    string("    Polynomial discretization for mixingPlane ribbons") +
                    string(" generator is not implemented")
                );
            }
            break;

            default:
            {
                FatalErrorInFunction
                    << "Bad type of mixing plane discretization: "
                    << discretizationMethodNames_[ribbonDiscretization_.method]
                    << nl << "Available types are: " << nl
                    << discretizationMethodNames_
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            InfoInFunction << setprecision(12) << nl
                << "masterBB: " << masterBB << nl
                << "slaveBB: " <<  slaveBB << nl
                << "globSpanBB: " <<  globSpanBB << nl
                << "initial profile values: " << profileBeforeClip << endl;
        }
    }

    // All done, let's transfer the profile to all processors
    Pstream::broadcast(profileBeforeClip);

    if(debug)
    {
        Pout << "MixingPlaneInterpolation::calcProfile(): profileBeforeClip values: "  << profileBeforeClip << endl;
    }

    // Allocate the return container
    tmp<pointField> tprofile(new pointField(profileBeforeClip.size()));
    pointField& profile = tprofile.ref();

#if 0
    // Martin to work here. HJ, 27/Jan/2011

    // Remove interpolationProfile_ points located outside of either
    // master/slave patch boundingBox,
    // with the exception of the first and last profile points

    // This editing of the profile points is not necessary.
    // The existing ribbon profile generators are already using the mixingPlane
    // patches bounding box as profile extrema, so there is no need to search
    // and remove any profile points located outside of globSpanBB.
    // The user-defined profile generator might introduce some profile points
    // outside of the bounding box, but we have to assume this is a
    // configuration chosen by the user.
    // So unless we have a real use case for this, let's leave this commented
    // MB, 28/Mar/2021

    label curIndex = 0;

    forAll (profileBeforeClip, pI)
    {
        if (globSpanBB.contains(profileBeforeClip[pI]))
        {
            // We keep that profile point
            profile[curIndex] = profileBeforeClip[pI];
            curIndex++;  // Next slot
        }
        else
        {
            if (debug)
            {
                InfoIn
                (
                    "MixingPlaneInterpolation"
                    "<MasterPatch, SlavePatch>::"
                    "removeNonOverlappedProfilePoints"
                )   << setprecision(12) << "   Removing point: "
                    << profileBeforeClip[pI] << endl;
            }
        }
    }

    profile.setSize(curIndex);

#else

    // Transfer all profile points, without exception
    forAll (profileBeforeClip, pI)
    {
        profile[pI] = profileBeforeClip[pI];
    }

#endif

    if (profile.size() < 2)
    {
        FatalErrorIn
        (
            "tmp<pointField> MixingPlaneInterpolation<MasterPatch, "
            "SlavePatch>::calcProfile() const"
        )   << "Lost all points in profile: " << profile
            << abort(FatalError);
    }

    if (debug)
    {
        InfoIn
        (
            "MixingPlaneInterpolation<MasterPatch, SlavePatch>::"
            "removeNonOverlappedProfilePoints"
        )   << "cleaned-up profile values: " << profile << endl;
    }

    return tprofile;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
const Foam::pointField&
MixingPlaneInterpolation<MasterPatch, SlavePatch>::interpolationProfile() const
{
    if (interpolationProfile_.size() == 0)
    {
        // Compute profile or read user-defined profile from file
        interpolationProfile_ = calcProfile();
    }

    return interpolationProfile_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
