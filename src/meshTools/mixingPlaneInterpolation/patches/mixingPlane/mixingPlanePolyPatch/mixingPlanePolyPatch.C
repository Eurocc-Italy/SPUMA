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

Author
    Martin Beaudoin, Hydro-Quebec, 2009. All rights reserved.

Contributor
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "mixingPlanePolyPatch.H"
#include "polyBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "ZoneIDs.H"
#include "Time.H"
#include "SubField.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, dictionary);

    scalar
    Foam::mixingPlanePolyPatch::tolPointDistanceFromRzero_
    (
        debug::floatOptimisationSwitch("mixingPlaneTolPointDistanceFromRzero", 1.0e-10)
    );

    registerOptSwitch
    (
        "mixingPlaneTolPointDistanceFromRzero",
        scalar,
        mixingPlanePolyPatch::tolPointDistanceFromRzero_
    );

    bool
    Foam::mixingPlanePolyPatch::writeRibbonsProfileToFile_
    (
        debug::infoSwitch("writeMixingPlaneRibbonsProfileToFiles", false)
    );

    registerInfoSwitch
    (
        "writeMixingPlaneRibbonsProfileToFiles",
        bool,
        mixingPlanePolyPatch::writeRibbonsProfileToFile_
    );

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::mixingPlanePolyPatch::active() const
{
    polyPatchID shadow(shadowName_, boundaryMesh());
    faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

    // For decomposition and reconstruction
    // If not runing in parallel and the patch is not local, this is a serial
    // operation on a piece of a parallel decomposition and is therefore
    // inactive.  HJ, 5/Sep/2016
#ifdef FOAM_EXTEND_PARALLEL
    if (!Pstream::parRun() && !localParallel())
#else
    if (!false && !localParallel())
#endif
    {
        return false;
    }

    return shadow.active() && zone.active();
}


void Foam::mixingPlanePolyPatch::calcZoneAddressing() const
{
    // Calculate patch-to-zone addressing
    if (zoneAddressingPtr_)
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcZoneAddressing() const")
            << "Patch to zone addressing already calculated"
            << abort(FatalError);
    }

    // Calculate patch-to-zone addressing
    zoneAddressingPtr_ = new labelList(size());
    labelList& zAddr = *zoneAddressingPtr_;
    const faceZone& myZone = zone();

    for (label i = 0; i < size(); i++)
    {
        zAddr[i] = myZone.whichFace(start() + i);
    }

    // Check zone addressing
    if (zAddr.size() > 0 && min(zAddr) < 0)
    {
        Info<< "myZone: " << myZone << nl
            << "my start and size: " << start() << " and " << size() << nl
            << "zAddr: " << zAddr << endl;

        FatalErrorIn("void mixingPlanePolyPatch::calcZoneAddressing() const")
            << "Problem with patch-to-zone addressing: some patch faces "
            << "not found in interpolation zone"
            << abort(FatalError);
    }
}


void Foam::mixingPlanePolyPatch::calcPatchToPatch() const
{
    // Create patch-to-patch interpolation
    if (patchToPatchPtr_)
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcPatchToPatch() const")
            << "Patch to patch interpolation already calculated"
            << abort(FatalError);
    }

    if (master())
    {
        if (debug)
        {
            Info<< "Creating mixingPlaneInterpolation for patch "
                << name() << " with shadow " << shadowName() << nl
                << "ribbonDiscretization method = "
                << mixingPlaneInterpolation::discretizationMethodNames_
                       [ribbonDiscretization_.method]
                << " sweepAxisType = "
                << mixingPlaneInterpolation::sweepAxisNames_[sweepAxisType_]
                << " " << sweepAxisType_
                << " stackAxisType = "
                << mixingPlaneInterpolation::stackAxisNames_[stackAxisType_]
                << " " << stackAxisType_
                << endl;
        }

        patchToPatchPtr_ =
            new mixingPlaneZoneInterpolation
            (
                zone()(),
                shadow().zone()(),
                csPtr_(),
                ribbonDiscretization_,
                sweepAxisType_,
                stackAxisType_,
                wfEvalMethod_,
                tolPointDistanceFromRzero_
            );

        // If using the info switch -info-switch mixingPlaneWriteRibbonsProfileToFile=1,
        // write the mixingPlane ribbon generator data as CSV and VTK files
        if (writeRibbonsProfileToFile_)
        {
            const polyMesh& mesh = boundaryMesh().mesh();

            patchToPatchPtr_->writeMixingPlaneRibbonGenerationData
            (
                mesh.time().path(),
                name(),
                shadowName()
            );
        }

        // Let's do some sanity checks right away.
        // For a valid mixingPlane interface, there should be no GGI/AMI weighting
        // factors with zero values, or no interface patch facets without any
        // neighbours.
        // This kind of problem can occur for various reasons: a user-defined
        // profile with missing ribbons, not covering the full size of the
        // mixingPlane patches in the stacking direction; some large orientation
        // differences between the patch ribbons normal and the normal of some
        // their neighbouring mixingPlane patches facets, which will activate a
        // GGI quick-reject test based on the angle cosine of those normals,
        // etc.
        // When this situation occurs, we basically get into a non-overlapped
        // GGI facets configuration, which will cause the solver to usually
        // crash with a somewhat cryptic FPE exception caused by a division by
        // zero for the evaluation of some of the fields across the mixingPlane
        // interface.
        // So let's catch this right away so we can inform the user with a
        // useful message, and stop the solver immediatly.
        //
        // Unless, of course, the master and shadow patches have a size of zero.

        if(this->size() > 0 && shadow().size() > 0)
        {

            bool foundLonelyFacets = false;

            const scalarListList *masterShadowProfileGGIWeightLists[4] =
            {
                &(patchToPatchPtr_->masterPatchToProfileWeights()),
                &(patchToPatchPtr_->masterProfileToPatchWeights()),
                &(patchToPatchPtr_->slavePatchToProfileWeights()),
                &(patchToPatchPtr_->slaveProfileToPatchWeights())
            };

            // Iterate over the 4 GGI weights lists, looking for zero neighbours or
            // zero-value GGI weights
            for(int i=0; i<4; i++)
            {
                const scalarListList& weightListList =
                    *(masterShadowProfileGGIWeightLists[i]);

                forAll(weightListList, j)
                {
                    if(weightListList[j].size() == 0)
                    {
                        // No neighbours
                        foundLonelyFacets = true;
                        break;
                    }

                    forAll(weightListList[j], k)
                    {
                        if(weightListList[j][k] == 0.0)
                        {
                            // zero-value GGI weight
                            foundLonelyFacets = true;
                            break;
                        }
                    }
                    if(foundLonelyFacets == true)
                    {
                        break;
                    }
                }
                if(foundLonelyFacets == true)
                {
                    break;
                }
            }

            if(foundLonelyFacets)
            {
                FatalErrorInFunction
                    << "MixingPlane interface error: ("
                    << name() << ", " << shadowName() << ")"
                    << endl
                    << "Detecting at least one mixingPlane ribbon or one "
                    << "mixingPlane patch facet without any neighbours."
                    << endl
                    << "This is a fatal error. Some possible solutions:"
                    << endl
                    << " - Review your mixingPlane interface definition"
                    << endl
                    << " - Adjust your mesh in the area of the mixingPlane patches"
                    << endl
                    << " - Adjust your user-defined profile to cover the full extent"
                    << " of the mixingPlane patches in the stacking direction"
                    << endl
                    << " - Adjust your runtime parameters like the GGIFeatureCosTol "
                    << "tolerance factor"
                    << endl
                    << abort(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcPatchToPatch() const")
            << "Attempting to create MixingPlaneInterpolation on a shadow"
            << abort(FatalError);
    }
}


void Foam::mixingPlanePolyPatch::calcReconFaceCellCentres() const
{
    if (reconFaceCellCentresPtr_)
    {
        FatalErrorIn
        (
            "void mixingPlanePolyPatch::calcReconFaceCellCentres() const"
        )   << "Reconstructed cell centres already calculated"
            << abort(FatalError);
    }

    // Create neighbouring face centres using interpolation
    if (master())
    {
        const label shadowID = shadowIndex();

        reconFaceCellCentresPtr_ =
            new vectorField
            (
                interpolate
                (
                    boundaryMesh()[shadowID].faceCellCentres()
                  - boundaryMesh()[shadowID].faceCentres()
                )
              + faceCentres()
            );
    }
    else
    {
        FatalErrorIn
        (
            "void mixingPlanePolyPatch::calcReconFaceCellCentres() const"
        )   << "Attempting to create reconFaceCellCentres on a shadow"
            << abort(FatalError);
    }
}


void Foam::mixingPlanePolyPatch::calcLocalParallel() const
{
    // Calculate patch-to-zone addressing
    if (localParallelPtr_)
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcLocalParallel() const")
            << "Local parallel switch already calculated"
            << abort(FatalError);
    }

    localParallelPtr_ = new bool(false);
    bool& emptyOrComplete = *localParallelPtr_;

    // If running in serial, all mixingPlanes are expanded to zone size.
    // This happens on decomposition and reconstruction where
    // size and shadow size may be zero, but zone size may not
    // HJ, 1/Jun/2011
#ifdef FOAM_EXTEND_PARALLEL
    if (!Pstream::parRun())
#else
    if (!false)
#endif
    {
        emptyOrComplete = false;
    }
    else
    {
        // Check that patch size is greater than the zone size.
        // This is an indication of the error where the face zone is not global
        // in a parallel run.  HJ, 9/Nov/2014
        if (size() > zone().size())
        {
            FatalErrorIn
            (
                "void mixingPlanePolyPatch::calcLocalParallel() const"
            )   << "Patch size is greater than zone size for GGI patch "
                << name() << ".  This is not allowerd: "
                << "the face zone must contain all patch faces and be "
                << "global in parallel runs"
                << abort(FatalError);
        }

        // Calculate localisation on master and shadow
        emptyOrComplete =
            (
                zone().size() == size()
             && shadow().zone().size() == shadow().size()
            )
         || (size() == 0 && shadow().size() == 0);

        reduce(emptyOrComplete, andOp<bool>());
    }

#ifdef FOAM_EXTEND_PARALLEL
    if (debug && Pstream::parRun())
#else
    if (debug && false)
#endif
    {
        Info<< "mixingPlane patch Master: " << name()
            << " Slave: " << shadowName() << " is ";

        if (emptyOrComplete)
        {
           Info<< "local parallel" << endl;
        }
        else
        {
            Info<< "split between multiple processors" << endl;
        }
    }
}


void Foam::mixingPlanePolyPatch::calcReceive() const
{
    // Note: all processors will execute calcReceive but only master will
    // hold the information.  Therefore, pointers on slave processors
    // will remain meaningless, but for purposes of consistency
    // (of the calc-call) they will be set to zero-sized array
    // HJ, 4/Jun/2011

    if (receiveAddrPtr_)
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcReceive() const")
            << "Receive addressing already calculated"
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "mixingPlanePolyPatch::calcReceive() const for patch "
            << index() << endl;
    }

#ifdef FOAM_EXTEND_PARALLEL
    if (!Pstream::parRun())
#else
    if (!false)
#endif
    {
        FatalErrorIn("void mixingPlanePolyPatch::calcReceive() const")
            << "Requested calculation of send-receive addressing for a "
            << "serial run.  This is not allowed"
            << abort(FatalError);
    }

    // Master will receive and store the maps
    if (Pstream::master())
    {
        receiveAddrPtr_ = new labelListList(Pstream::nProcs());
        labelListList& rAddr = *receiveAddrPtr_;

        // Insert master
        rAddr[0] = zoneAddressing();

        for (label procI = 1; procI < Pstream::nProcs(); procI++)
        {
            // Note: must use normal comms because the size of the
            // communicated lists is unknown on the receiving side
            // HJ, 4/Jun/2011

            // Opt: reconsider mode of communication
            IPstream ip(Pstream::commsTypes::scheduled, procI);

            rAddr[procI] = labelList(ip);
        }
    }
    else
    {
        // Create dummy pointers: only master processor stores maps
        receiveAddrPtr_ = new labelListList();

        // Send information to master
        const labelList& za = zoneAddressing();

        // Note: must use normal comms because the size of the
        // communicated lists is unknown on the receiving side
        // HJ, 4/Jun/2011

        // Opt: reconsider mode of communication
        OPstream op(Pstream::commsTypes::scheduled, Pstream::masterNo());

        // Send local and remote addressing to master
        op << za;
    }
}


const Foam::labelListList& Foam::mixingPlanePolyPatch::receiveAddr() const
{
    if (!receiveAddrPtr_)
    {
        calcReceive();
    }

    return *receiveAddrPtr_;
}


void Foam::mixingPlanePolyPatch::clearLocalGeom() const
{
    deleteDemandDrivenData(reconFaceCellCentresPtr_);
}


void Foam::mixingPlanePolyPatch::clearOut() const
{
    clearLocalGeom();

    shadowIndex_ = -1;
    zoneIndex_ = -1;

    // Receive maps do not depend on the local position
    deleteDemandDrivenData(receiveAddrPtr_);

    deleteDemandDrivenData(zoneAddressingPtr_);
    deleteDemandDrivenData(patchToPatchPtr_);
    deleteDemandDrivenData(localParallelPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, transform),
    neighbPatchName_(fileName::null),
    shadowName_(fileName::null),
    zoneName_("initializeMe"),
    neighbPatchID_(-1),
    coupleGroup_(),

    csPtr_
    (
        new fe::coordinateSystem
        (
            "mixingCS",
            vector::zero,
            vector(0, 0, 1),
            vector(1, 0, 0)
        )
    ),
    ribbonDiscretization_(),
    sweepAxisType_(mixingPlaneInterpolation::SWEEP_UNKNOWN),
    stackAxisType_(mixingPlaneInterpolation::STACK_UNKNOWN),
    wfEvalMethod_(mixingPlaneInterpolation::WFEVALMETHOD_UNKNOWN),
    userProfileFile_(fileName::null),
    writeAveragedPatchField_(false),
    shadowIndex_(-1),
    zoneIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    receiveAddrPtr_(nullptr)
{
}

Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    neighbPatchName_(dict.getOrDefault("neighbourPatch", word::null)),
    shadowName_(neighbPatchName_),
    zoneName_(dict.lookup("zone")),
    neighbPatchID_(-1),
    coupleGroup_(dict),
    csPtr_
    (
        new fe::coordinateSystem
        (
            "mixingCS",
            vector::zero,
            vector(0, 0, 1),
            vector(1, 0, 0)
        )
    ),
    ribbonDiscretization_(),
    sweepAxisType_(mixingPlaneInterpolation::SWEEP_UNKNOWN),
    stackAxisType_(mixingPlaneInterpolation::STACK_UNKNOWN),
    wfEvalMethod_(mixingPlaneInterpolation::WFEVALMETHOD_UNKNOWN),
    userProfileFile_(fileName::null),
    writeAveragedPatchField_(dict.getOrDefault("writeAveragedPatchField", false)),
    shadowIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    receiveAddrPtr_(nullptr)
{
    // Sanity checks
    if (neighbPatchName_ == word::null && !coupleGroup_.valid())
    {
        Info << "DEBUG: dict: " << dict << endl;

        FatalIOErrorInFunction(dict)
            << "No \"neighbourPatch\" provided." << endl
            << "Is your mesh uptodate with split ggis?" << endl
            << "Run foamUpgradeCyclics to convert mesh and fields"
            << " to split ggis." << exit(FatalIOError);
    }

    if (neighbPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch name " << neighbPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // When constructing from dictionary, only master side information will be
    // read and used.  This requires special check, because polyBoundaryMesh
    // is being filled at this point.  See fix in findPatchID.
    // HJ, and MB, 28/Jan/2011

    // Check if shadow exists.  If so, we are on the slave side
    Foam::label shadowPatchId = boundaryMesh().findPatchID(shadowName_);

    if (shadowPatchId == -1)
    {
        // Master side, read additional data

        csPtr_ =
            fe::coordinateSystem::New
            (
                "mixingCS",
                dict.subDict("coordinateSystem")
            );

        const dictionary& ribbonPatchSubDict = dict.subDict("ribbonPatch");

        sweepAxisType_ =
            mixingPlaneInterpolation::sweepAxisNames_.read
            (
                ribbonPatchSubDict.lookup("sweepAxis")
            );

        stackAxisType_ =
            mixingPlaneInterpolation::stackAxisNames_.read
            (
                ribbonPatchSubDict.lookup("stackAxis")
            );

        wfEvalMethod_ =
            mixingPlaneInterpolation::wfEvalMethodNames_
            [
                dict.lookupOrDefault<word>("weightFactorsEvalMethod", "ggi")
            ];

        ribbonDiscretization_.method =
            mixingPlaneInterpolation::discretizationMethodNames_.read
            (
                ribbonPatchSubDict.lookup("discretization")
            );

        if (ribbonDiscretization_.method == mixingPlaneInterpolation::FOLLOW_PATCH_DISCRETIZATION)
        {
            // Memorize followPatchDiscretization params
            ribbonDiscretization_.methodParamsDict =
                ribbonPatchSubDict.subDict("followPatchDiscretizationParams");

            // Sanity check: verify the optional histogram's bin size scale
            // factor is greater than zero
            scalar histoBinSizeScalingFactor =
                ribbonDiscretization_.methodParamsDict.lookupOrDefault<scalar>
                (
                    "histoBinSizeScaleFactor",
                    1.0
                );

            if (histoBinSizeScalingFactor <= 0.0)
            {
                FatalErrorInFunction
                    << "mixingPlane interface : " << name
                    << "/" << shadowName()  << endl
                    << "The scaling factor specified by the dictionary entry "
                    << "'histoBinSizeScaleFactor' cannot be zero or negative."
                    << abort(FatalError);
            }


        }

        else if (ribbonDiscretization_.method == mixingPlaneInterpolation::LINEAR_ALONG_STACK_AXIS)
        {
            // Memorize linear params
            ribbonDiscretization_.methodParamsDict =
                ribbonPatchSubDict.subDict("linearAlongStackAxisParams");

            // Validation
            label nRibbons =
                ribbonDiscretization_.methodParamsDict.lookupOrDefault<label>("nRibbons", 0);

            if(nRibbons <= 0)
            {
                FatalErrorInFunction
                    << "Error in file constant/polyMesh/boundary option:" << nl
                    << "    MixingPlane patches: "  <<  zoneName_ << "/" << shadowName_ << nl
                    << "    MixingPlane discretization method: linearAlongStackAxis" << nl
                    << "    Parameter value 'nRibbons' must be greater than 0"
                    << abort(FatalError);
            }
        }

        else if (ribbonDiscretization_.method == mixingPlaneInterpolation::POLYNOMIAL)
        {
            // Memorize polynomial params
            ribbonDiscretization_.methodParamsDict =
                ribbonPatchSubDict.subDict("polynomialParams");
        }

        else if (ribbonDiscretization_.method == mixingPlaneInterpolation::USER_DEFINED)
        {
            // Profile coordinates are provided through a CSV file.
            // Memorize CSV file params
            ribbonDiscretization_.methodParamsDict =
                ribbonPatchSubDict.subDict("userDefinedProfileParams");

            dictionary& dict = ribbonDiscretization_.methodParamsDict;

            // Make sure the CSV file is read from the constant directory
            fileName csvFile(dict.lookup("file"));
            dict.add
            (
                "file",
                boundaryMesh().mesh().time().constantPath()/csvFile.name(),
                true  // merge or replace
            );

            // Check for non-implemented features

            // Check for interpolation of profile
            label nbrRequestedRibbons = dict.lookupOrDefault("nRibbons", 0);

            if(nbrRequestedRibbons > 0)
            {
                notImplemented
                (
                    string(FUNCTION_NAME) + nl +
                    string("    Error in mixingPlane patch definition: ") +
                    zoneName_ + string("/") + shadowName_ + nl +
                    string("    In section: userDefinedProfileParams:") + nl +
                    string("    Parameter: nRibbons ") + std::to_string(nbrRequestedRibbons) + nl +
                    string("    Interpolation or sampling of user-defined mixingPlane") +
                    string(" profile is not implemented")
                );
            }


        }
    }
}


Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    neighbPatchName_(pp.neighbPatchName_),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    neighbPatchID_(pp.neighbPatchID_),
    coupleGroup_(pp.coupleGroup_),
    csPtr_(pp.csPtr_->clone()),
    ribbonDiscretization_(pp.ribbonDiscretization_),
    sweepAxisType_(pp.sweepAxisType_),
    stackAxisType_(pp.stackAxisType_),
    wfEvalMethod_(pp.wfEvalMethod_),
    userProfileFile_(pp.userProfileFile_),
    writeAveragedPatchField_(pp.writeAveragedPatchField_),
    shadowIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    receiveAddrPtr_(nullptr)
{
}

//- Copy construct, resetting nbrPatchID and faceCells
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    label nbrPatchID,
    const labelList& faceCells
)
:
    coupledPolyPatch(pp, faceCells),
    neighbPatchName_(pp.neighbPatchName_),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    neighbPatchID_(nbrPatchID), // Or pp.neighbPatchID_ ??
    coupleGroup_(pp.coupleGroup_),
    csPtr_(pp.csPtr_->clone()),
    ribbonDiscretization_(pp.ribbonDiscretization_),
    sweepAxisType_(pp.sweepAxisType_),
    stackAxisType_(pp.stackAxisType_),
    wfEvalMethod_(pp.wfEvalMethod_),
    userProfileFile_(pp.userProfileFile_),
    writeAveragedPatchField_(pp.writeAveragedPatchField_),
    shadowIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    receiveAddrPtr_(nullptr)

{
    InfoInFunction << "Warning: Need to verify if we need this constructor."  << endl;
}

// Construct as copy, resetting the face list and boundary mesh data
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& neighbPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    neighbPatchName_(neighbPatchName),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    neighbPatchID_(pp.neighbPatchID_),
    coupleGroup_(pp.coupleGroup_),
    csPtr_(pp.csPtr_->clone()),
    ribbonDiscretization_(pp.ribbonDiscretization_),
    sweepAxisType_(pp.sweepAxisType_),
    stackAxisType_(pp.stackAxisType_),
    wfEvalMethod_(pp.wfEvalMethod_),
    userProfileFile_(pp.userProfileFile_),
    writeAveragedPatchField_(pp.writeAveragedPatchField_),
    shadowIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    receiveAddrPtr_(nullptr)
{
}


Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    neighbPatchName_(pp.neighbPatchName_),
    shadowName_(pp.shadowName_),
    zoneName_(pp.zoneName_),
    neighbPatchID_(pp.neighbPatchID_),
    coupleGroup_(pp.coupleGroup_),
    csPtr_(pp.csPtr_->clone()),
    ribbonDiscretization_(pp.ribbonDiscretization_),
    sweepAxisType_(pp.sweepAxisType_),
    stackAxisType_(pp.stackAxisType_),
    wfEvalMethod_(pp.wfEvalMethod_),
    userProfileFile_(pp.userProfileFile_),
    writeAveragedPatchField_(pp.writeAveragedPatchField_),
    shadowIndex_(-1),
    patchToPatchPtr_(nullptr),
    zoneAddressingPtr_(nullptr),
    reconFaceCellCentresPtr_(nullptr),
    localParallelPtr_(nullptr),
    receiveAddrPtr_(nullptr)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingPlanePolyPatch::~mixingPlanePolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::mixingPlanePolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1 && shadowName_ != Foam::word::null)
    {
        // Grab shadow patch index
        polyPatchID shadow(shadowName_, boundaryMesh());

        if (!shadow.active())
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadowName_
                << " not found.  Please check your MixingPlane definition.  "
                << "This may be fine at mesh generation stage."
                << endl;
        }

        shadowIndex_ = shadow.index();

        // Check the other side is a mixingPlane
        if (!isA<mixingPlanePolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "Shadow of mixingPlane patch " << name()
                << " named " << shadowName() << " is not a mixingPlane.  "
                << "Type: " << boundaryMesh()[shadowIndex_].type() << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }

        // Check for mixingPlane onto self
        if (index() == shadowIndex_)
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "mixingPlane patch " << name()
                << " created as its own shadow"
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}


const Foam::fe::coordinateSystem& Foam::mixingPlanePolyPatch::cs() const
{
    if (master())
    {
        return csPtr_();
    }
    else
    {
        return shadow().cs();
    }
}


Foam::label Foam::mixingPlanePolyPatch::zoneIndex() const
{
    if (zoneIndex_ == -1 && zoneName_ != Foam::word::null)
    {
        // Grab zone patch index
        faceZoneID zone(zoneName_, boundaryMesh().mesh().faceZones());

        if (!zone.active())
        {
            FatalErrorIn("label mixingPlanePolyPatch::zoneIndex() const")
                << "Face zone name " << zoneName_
                << " for mixingPlane patch " << name()
                << " not found.  "
                << "Please check mixingPlane interface definition."
                << abort(FatalError);
        }

        zoneIndex_ = zone.index();
    }

    return zoneIndex_;
}


const Foam::mixingPlanePolyPatch& Foam::mixingPlanePolyPatch::shadow() const
{
    return refCast<const mixingPlanePolyPatch>(boundaryMesh()[shadowIndex()]);
}


const Foam::faceZone& Foam::mixingPlanePolyPatch::zone() const
{
    return boundaryMesh().mesh().faceZones()[zoneIndex()];
}


const Foam::labelList& Foam::mixingPlanePolyPatch::zoneAddressing() const
{
    if (!zoneAddressingPtr_)
    {
        calcZoneAddressing();
    }

    return *zoneAddressingPtr_;
}


bool Foam::mixingPlanePolyPatch::localParallel() const
{
    // Calculate patch-to-zone addressing
    if (!localParallelPtr_)
    {
        calcLocalParallel();
    }

    return *localParallelPtr_;
}


const Foam::mixingPlaneZoneInterpolation&
Foam::mixingPlanePolyPatch::patchToPatch() const
{
    if (master())
    {
        if (!patchToPatchPtr_)
        {
            Info<< "Initializing the mixingPlane interpolator between "
                << "master/shadow patches: "
                << name() << "/" << shadowName()
                << endl;


            calcPatchToPatch();
        }

        return *patchToPatchPtr_;
    }
    else
    {
        return shadow().patchToPatch();
    }
}


Foam::label Foam::mixingPlanePolyPatch::comm() const
{
    return boundaryMesh().mesh().comm();
}


int Foam::mixingPlanePolyPatch::tag() const
{
    return Pstream::msgType();
}


const Foam::vectorField&
Foam::mixingPlanePolyPatch::reconFaceCellCentres() const
{
    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }

    return *reconFaceCellCentresPtr_;
}

void Foam::mixingPlanePolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    // Communication is allowed either before or after processor
    // patch comms.  HJ, 11/Jul/2011
    if (active())
    {
        // Calculate transforms for correct mixingPlane cut
        calcTransforms();

        // Force zone addressing and remote zone addressing
        // (uses mixingPlane interpolator)
        zoneAddressing();

        // Force local parallel
#ifdef FOAM_EXTEND_PARALLEL
        if (Pstream::parRun() && !localParallel())
#else
        if (false && !localParallel())
#endif
        {
            // Calculate send addressing
            receiveAddr();
        }

        // Note: Only master calculates recon; slave uses master interpolation
        if (master())
        {
            reconFaceCellCentres();
        }
    }

    polyPatch::initGeometry(pBufs);
}

const Foam::word& Foam::mixingPlanePolyPatch::neighbPatchName() const
{
    if (neighbPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        neighbPatchName_ = boundaryMesh()[patchID].name();
    }
    return neighbPatchName_;
}

Foam::label Foam::mixingPlanePolyPatch::neighbPatchID() const
{
    if (neighbPatchID_ == -1)
    {
        neighbPatchID_ = this->boundaryMesh().findPatchID(neighbPatchName());

        if (neighbPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << neighbPatchName()
                << endl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a mixingPlane interface
        const mixingPlanePolyPatch& nbrPatch = refCast<const mixingPlanePolyPatch>
        (
            this->boundaryMesh()[neighbPatchID_]
        );

        if (nbrPatch.neighbPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << neighbPatchName()
                << endl << " but that in return specifies "
                << nbrPatch.neighbPatchName()
                << endl;
        }
    }
    return neighbPatchID_;
}

void Foam::mixingPlanePolyPatch::transformPosition(pointField& l) const
{
    WarningInFunction << "Need to verify if we need this"  << endl;
}

void Foam::mixingPlanePolyPatch::transformPosition
(
    point& l,
    const label facei
) const
{
    WarningInFunction << "Need to verify if we need this"  << endl;
}


void Foam::mixingPlanePolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{
}

void Foam::mixingPlanePolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    polyPatch::calcGeometry(pBufs);

    // Note: Calculation of transforms must be forced before the
    // reconFaceCellCentres in order to correctly set the transformation
    // in the interpolation routines.
    // HJ, 3/Jul/2009
}


void Foam::mixingPlanePolyPatch::initMovePoints(PstreamBuffers& pBufs, const pointField& p)
{
    clearLocalGeom();

    // Calculate transforms on mesh motion?
    calcTransforms();

    if (master())
    {
        shadow().clearLocalGeom();
        shadow().calcTransforms();
    }

    // Update interpolation for new relative position of mixingPlane interfaces
    if (patchToPatchPtr_)
    {
        patchToPatchPtr_->movePoints();
    }

    // Recalculate send and receive maps
    if (active())
    {
        // Force zone addressing first
        zoneAddressing();

#ifdef FOAM_EXTEND_PARALLEL
        if (Pstream::parRun() && !localParallel())
#else
        if (false && !localParallel())
#endif
        {
            receiveAddr();
        }
    }

    if (active() && master())
    {
        reconFaceCellCentres();
    }

    polyPatch::initMovePoints(pBufs, p);
}


void Foam::mixingPlanePolyPatch::movePoints(PstreamBuffers& pBufs, const pointField& p)
{
    polyPatch::movePoints(pBufs, p);
}


void Foam::mixingPlanePolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::mixingPlanePolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    clearOut();
}


void Foam::mixingPlanePolyPatch::calcTransforms() const
{
}


void Foam::mixingPlanePolyPatch::initOrder(PstreamBuffers& pBufs, const primitivePatch& pp) const
{}


bool Foam::mixingPlanePolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size(), -1);
    rotation.setSize(pp.size(), 0);

    // Nothing changes
    return false;
}


void Foam::mixingPlanePolyPatch::syncOrder() const
{}

bool Foam::mixingPlanePolyPatch::checkGeometry
(
    const bool report
) const
{
    return false;
}

bool Foam::mixingPlanePolyPatch::checkTopology
(
    const bool report
) const
{
    // Return value: return false if no problem detected
    bool foundProblem = false;

    // Checking for mixingPlane patch facets intersecting the rotation axis
    // r=0 without an explicit discretization point.
    // Will generate a modifyMeshDict for triangulating this patch facet.
    if(debug)
    {
        Info << "    Checking mixingPlane patch : " << name()
            << " for facet intersecting the r=0 axis.";
    }

    bool foundDiscretizedRZero = false;

    const pointField patchLocalCoords(cs().localPosition(localPoints()));

    scalar tolDistanceFromRzero = fabs(tolPointDistanceFromRzero_);

    // Check facets discretization first
    forAll (patchLocalCoords, i)
    {
        scalar r = fabs(patchLocalCoords[i][0]);

        if(r <= tolDistanceFromRzero)
        {
            // All done, we have at least one discretization point from this
            // patch ar r=0.
            foundDiscretizedRZero = true;
            break;
        }
    }

    if (!foundDiscretizedRZero)
    {
        // No discretization point were found at r=0. Let's find facets possibly
        // intersecting the rotation axis e3.

        bool foundIntersect = false;
        point pointAtRzero;
        label localFaceIndex = 0;

        // Iterate through patch faces
        const faceList& patchFaces = localFaces();

        forAll (patchFaces, i)
        {
            // Search intersection with rotation axis e3
            pointHit inter = patchFaces[i].ray
                (
                    cs().origin(),
                    cs().e3(),
                    localPoints()
                );

            if (inter.hit())
            {
                localFaceIndex = i;
                pointAtRzero = inter.hitPoint();
                foundIntersect = true;
                break;
            }
        }

        if (foundIntersect)
        {
            foundProblem = true;

            label newSplitFaceId = localFaceIndex + start();

            if(report)
            {
                Info << "\n*** Not OK. MixingPlane patch "
                    << name()
                    << " requires a discretization point at r=0: "
                    << pointAtRzero
                    << " for face index: " << newSplitFaceId
                    << "\n*** Writing a new modifyMeshDict entry for fixing"
                    << " this issue" << endl;


                const polyMesh& mesh = boundaryMesh().mesh();

                IOdictionary dict
                (
                    IOobject
                    (
                        "modifyMeshDict",
                        mesh.time().system(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    )
                );

                // Check required entries, add empty entries if none found

                std::list<word> requiredEntries = {
                    "pointsToMove",
                    "edgesToSplit",
                    "facesToTriangulate",
                    "faceLabelsToTriangulate",
                    "edgesToCollapse",
                    "cellsToSplit"
                };

                typedef List<word>  emptyEntry;

                for(auto e : requiredEntries)
                {
                    if (!dict.found(e))
                    {
                        dict.add(e, emptyEntry());
                    }
                }

                // Grab current list
                List<Tuple2<label, point> > faceLabelsToTriangulate
                (
                    dict.lookup("faceLabelsToTriangulate")
                );

                // Add new face/point combo, avoiding duplicate entries
                bool foundDuplicate = false;

                forAll(faceLabelsToTriangulate, i)
                {
                    if(faceLabelsToTriangulate[i].first() == newSplitFaceId)
                    {
                        foundDuplicate = true;
                        break;
                    }
                }

                if(!foundDuplicate)
                {
                    Tuple2<label, point> newFacePoint
                    (
                        newSplitFaceId,
                        pointAtRzero
                    );

                    faceLabelsToTriangulate.append(newFacePoint);
                }

                // Update dictionary entry
                dict.add("faceLabelsToTriangulate", faceLabelsToTriangulate, true);

                // Write updated modifyMeshDict
                dict.regIOobject::write();
            }
        }
    }

    return foundProblem;
}


void Foam::mixingPlanePolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("shadowPatch") << shadowName_ << token::END_STATEMENT << nl;
    os.writeKeyword("neighbourPatch") << neighbPatchName_ << token::END_STATEMENT << nl;
    os.writeKeyword("zone") << zoneName_ << token::END_STATEMENT << nl;

    // Note: only master writes the data
    if (master() || shadowIndex_ == -1)
    {
        // Write coordinate system dictionary.  Check by hand.  HJ, 26/Jan/2011
        os.writeKeyword("coordinateSystem");
        csPtr_().writeDict(os, true);

        dictionary ribbonPatchDict("ribbonPatch");

        ribbonPatchDict.add
        (
            "sweepAxis",
            mixingPlaneInterpolation::sweepAxisNames_[sweepAxisType_]
        );

        ribbonPatchDict.add
        (
            "stackAxis",
            mixingPlaneInterpolation::stackAxisNames_[stackAxisType_]
        );

        ribbonPatchDict.add
        (
            "discretization",
            mixingPlaneInterpolation::discretizationMethodNames_[ribbonDiscretization_.method]
        );

        switch (ribbonDiscretization_.method)
        {
            case MixingPlaneInterpolationName::FOLLOW_PATCH_DISCRETIZATION:
            {
                ribbonPatchDict.add
                (
                    "followPatchDiscretizationParams",
                    ribbonDiscretization_.methodParamsDict
                );
            }
            break;

            case MixingPlaneInterpolationName::LINEAR_ALONG_STACK_AXIS:
            {
                ribbonPatchDict.add
                (
                    "linearAlongStackAxisParams",
                    ribbonDiscretization_.methodParamsDict
                );
            }
            break;

            case MixingPlaneInterpolationName::POLYNOMIAL:
            {
                ribbonPatchDict.add
                (
                    "polynomialParams",
                    ribbonDiscretization_.methodParamsDict
                );
            }
            break;

            case MixingPlaneInterpolationName::USER_DEFINED:
            {
                ribbonPatchDict.add
                (
                    "userDefinedProfileParams",
                    ribbonDiscretization_.methodParamsDict
                );
            }
            break;

            // Deprecated but still supported
            case MixingPlaneInterpolationName::MASTER_PATCH:
            case MixingPlaneInterpolationName::SLAVE_PATCH:
            case MixingPlaneInterpolationName::BOTH_PATCHES:
            default:
            break;
        }

        os.writeKeyword("ribbonPatch")
            << ribbonPatchDict << nl;

        if (userProfileFile_ != fileName::null)
        {
            os.writeKeyword("userProfileFile") << userProfileFile_
                << token::END_STATEMENT << nl;
        }

        os.writeKeyword("weightFactorsEvalMethod")
            << mixingPlaneInterpolation::wfEvalMethodNames_[wfEvalMethod_]
            << token::END_STATEMENT << nl;

        os.writeKeyword("writeAveragedPatchField") << writeAveragedPatchField_
                << token::END_STATEMENT << nl;


    }
}


// ************************************************************************* //
