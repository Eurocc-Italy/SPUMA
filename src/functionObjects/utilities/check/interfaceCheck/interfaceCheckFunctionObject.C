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
    Ettore Fadiga, (2025)

\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\

    Copyright (C) 2023 Martin Beaudoin
    Copyright (C) 2025 CINECA (Ettore Fadiga)

-------------------------------------------------------------------------------

Class
    Foam::interfaceCheckFunctionObject

Description
    FunctionObject checking the fluxes across a pair of ggi, AMI or mixingPlane
    patches

Author
    Martin Beaudoin, 2023. All right reserved
    Ettore Fadiga, 2025. All right reserved

\*---------------------------------------------------------------------------*/

#include "interfaceCheckFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
// #include "mixingPlaneFvsPatchFields.H"
#include "overlapAMIFvsPatchFields.H"
// #include "ggiFvsPatchFields.H"
#include "cyclicAMIFvsPatchFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

    defineTypeNameAndDebug(interfaceCheckFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        interfaceCheckFunctionObject,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceCheckFunctionObject::interfaceCheckFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    phiName_(dict.lookup("phi"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info << "Creating interfaceCheck check function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceCheckFunctionObject::start()
{
    return true;
}


bool Foam::functionObjects::interfaceCheckFunctionObject::execute()
{
    const objectRegistry& mesh =
        time_.lookupObject<objectRegistry>(regionName_);

    if (mesh.foundObject<surfaceScalarField>(phiName_))
    {
        Info<< "Checking flux '" << phiName_ << "' balance across non-conformal boundary patches:" << endl;
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        boolList visited(phi.boundaryField().size(), false);

        label masterPatchI(-1);
        label neighbourPatchI(-1);
        word masterPatchName;
        word neighbourPatchName;
        word interfaceTypeName;

        forAll (phi.boundaryField(), patchI)
        {
            if (visited[patchI])
            {
                continue;
            }
            /*else if (isA<ggiFvsPatchScalarField>(phi.boundaryField()[patchI]))
            {
                const ggiPolyPatch& ggiPatch =
                    refCast<const ggiPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                masterPatchI = patchI;
                neighbourPatchI = ggiPatch.neighbPatchID();
                masterPatchName = ggiPatch.name();
                neighbourPatchName = ggiPatch.neighbPatchName();
                interfaceTypeName = ggiPatch.typeName;

                visited[masterPatchI] = true;
                visited[neighbourPatchI] = true;
            }*/
            /*else if (isA<mixingPlaneFvsPatchScalarField>(phi.boundaryField()[patchI]))
            {
                const mixingPlanePolyPatch& mixingPlanePatch =
                    refCast<const mixingPlanePolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                masterPatchI = patchI;
                neighbourPatchI = mixingPlanePatch.neighbPatchID();
                masterPatchName = mixingPlanePatch.name();
                neighbourPatchName = mixingPlanePatch.neighbPatchName();
                interfaceTypeName = mixingPlanePatch.typeName;

                visited[masterPatchI] = true;
                visited[neighbourPatchI] = true;
            }*/
            else if (isA<overlapAMIFvsPatchScalarField>(phi.boundaryField()[patchI]))
            {
                const overlapAMIPolyPatch& overlapAMIPatch =
                    refCast<const overlapAMIPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                masterPatchI = patchI;
                neighbourPatchI = overlapAMIPatch.neighbPatchID();
                masterPatchName = overlapAMIPatch.name();
                neighbourPatchName = overlapAMIPatch.neighbPatchName();
                interfaceTypeName = overlapAMIPatch.typeName;

                visited[masterPatchI] = true;
                visited[neighbourPatchI] = true;
            }

            else if (isA<cyclicAMIFvsPatchScalarField>(phi.boundaryField()[patchI]))
            {
                const cyclicAMIPolyPatch& cyclicAMIPatch =
                    refCast<const cyclicAMIPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                masterPatchI = patchI;
                neighbourPatchI = cyclicAMIPatch.neighbPolyPatchID();
                masterPatchName = cyclicAMIPatch.name();
                neighbourPatchName = cyclicAMIPatch.neighbPatchName();
                interfaceTypeName = cyclicAMIPatch.typeName;

                visited[masterPatchI] = true;
                visited[neighbourPatchI] = true;
            }

            if(masterPatchI != -1)
            {
                // Calculate flux and area from master and neighbour patches
                const scalar masterPatchArea =
                    gSum(phi.mesh().magSf().boundaryField()[masterPatchI]);

                const scalar masterPatchFlux = gSum(phi.boundaryField()[masterPatchI]);
                const scalar masterPatchFluxMag =
                    gSum(mag(phi.boundaryField()[masterPatchI]));

                // Calculate neighbour flux and area
                const scalar neighbourPatchArea =
                    gSum(phi.mesh().magSf().boundaryField()[neighbourPatchI]);

                const scalar neighbourPatchFlux =
                    gSum(phi.boundaryField()[neighbourPatchI]);
                const scalar neighbourPatchFluxMag =
                    gSum(mag(phi.boundaryField()[neighbourPatchI]));

                // Report
                Info << "Interface pair (" << masterPatchName << ", "
                    << neighbourPatchName << ")" << " : "
                    << "Interface type: " << interfaceTypeName << nl
                    << " Area: " << masterPatchArea << " " << neighbourPatchArea
                    << " Diff = " << masterPatchArea - neighbourPatchArea << " or "
                    << mag(masterPatchArea - neighbourPatchArea)/
                    (Foam::max(masterPatchArea, neighbourPatchArea) + SMALL)*100
                    << " %" << nl
                    << " Flux: " << masterPatchFluxMag << " " << neighbourPatchFluxMag
                    << " Diff = " << masterPatchFlux + neighbourPatchFlux << " or "
                    << mag(masterPatchFlux + neighbourPatchFlux)/
                    (Foam::max(masterPatchFluxMag, neighbourPatchFluxMag) + SMALL)*100
                    << " %" << endl;

                // Report scaled flux for mixingPlane interfaces
                if(interfaceTypeName == "mixingPlane" || interfaceTypeName == "overlapAMI" || interfaceTypeName == "cyclicAMI")
                {
                    scalar neighbourPatchScaleFactor_ = masterPatchArea/neighbourPatchArea;
                    scalar neighbourPatchFluxScaled  = neighbourPatchScaleFactor_ * neighbourPatchFlux;
                    scalar neighbourPatchFluxMagScaled = neighbourPatchScaleFactor_ * neighbourPatchFluxMag;

                    Info << " Scaled flux: "<< masterPatchFluxMag << " " << neighbourPatchFluxMagScaled
                        << " Diff = " << masterPatchFlux + neighbourPatchFluxScaled << " or "
                        << mag(masterPatchFlux + neighbourPatchFluxScaled)/
                        (Foam::max(masterPatchFluxMag, neighbourPatchFluxMagScaled) + SMALL)*100
                        << " %" << endl;
                }
            }

            // Reset patch numbers in case of multiples GGI or AMI interfaces
            masterPatchI = -1;
            neighbourPatchI = -1;
        }

        Info << endl;

        return true;
    }
    else
    {
        InfoIn("bool interfaceCheckFunctionObject::execute()")
            << "Cannot find flux field " << phiName_
            << endl;

        return false;
    }
}


bool Foam::functionObjects::interfaceCheckFunctionObject::timeSet()
{
    return true;
}


bool Foam::functionObjects::interfaceCheckFunctionObject::read(const dictionary& dict)
{
    return false;
}

// ************************************************************************* //
