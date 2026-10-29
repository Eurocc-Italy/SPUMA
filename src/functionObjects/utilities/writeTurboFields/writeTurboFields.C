/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
    Copyright (C) 2025 CINECA (Ettore Fadiga)
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

#include "Time.H"
#include "polyMesh.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"
#include "writeTurboFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeTurboFields, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeTurboFields,
        dictionary
    );
}
}

const Foam::objectRegistry& setObjRegistry
(
    const Foam::Time& runTime,
    const Foam::dictionary& dict
)
{
    const Foam::word regionName =
        dict.getOrDefault("region", Foam::polyMesh::defaultRegion);

    if (regionName == "__TIME__")
    {
        return runTime;
    }

    return runTime.lookupObject<Foam::objectRegistry>(regionName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeTurboFields::writeTurboFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name,runTime,dict),
    obr_(setObjRegistry(runTime, dict)),
    objectNames_(),
    pName_(""),
    UName_(""),
    TName_(""),
    MRF_(mesh_)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeTurboFields::read(const dictionary& dict)
{
    pName_ = dict.lookupOrDefault<word>("pName", "p");
    UName_ = dict.lookupOrDefault<word>("UName", "U");
    TName_ = dict.lookupOrDefault<word>("TName", "T");

    if (fvMeshFunctionObject::read(dict))
    {
        if (dict.found("field"))
        {
            objectNames_.resize(1);
            dict.readEntry("field", objectNames_.first());
        }
        else if (dict.found("fields"))
        {
            dict.readEntry("fields", objectNames_);
        }
        else
        {
            dict.readEntry("objects", objectNames_);
        }

        axialDirection_ = normalised(dict.lookupOrDefault<vector>("axialDir", vector(1.0,0.0,0.0)));
        origin_ = normalised(dict.lookupOrDefault<vector>("origin", vector(0.0,0.0,0.0)));
        omega_ = dict.lookupOrDefault<scalar>("omega", 0.0);
    }

    return true;
}


bool Foam::functionObjects::writeTurboFields::execute()
{
    return true;
}


bool Foam::functionObjects::writeTurboFields::write()
{
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);
    const volScalarField& T = obr_.lookupObject<volScalarField>(TName_);

    volVectorField URel("URel",U);
    MRF_.makeRelative(URel);
    volVectorField Utr("Utr", U - URel);

    const volScalarField::Boundary&  Tp = T.boundaryField();
    const basicThermo& thermo = basicThermo::lookupThermo(Tp[0]);
    const volScalarField&  h = thermo.he();
    tmp<volScalarField> rhot = thermo.rho();
    const volScalarField& rho = rhot.cref();

    // Finding radial and tangential direction for each patch face center
    const symmTensor planeDir
    (
        1 - axialDirection_.x(), 
        0.0,           
        0.0,
        1 - axialDirection_.y(), 
        0.0,
        1 - axialDirection_.z()
    );

    volScalarField UMer
    (
        IOobject
        (
            "UMer",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        Zero,
        dimVelocity,
        "calculated"
    );

    volScalarField URelMer
    (
        IOobject
        (
            "URelMer",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        Zero,
        dimVelocity, 
        "calculated"
    );

    volVectorField Ucyl
    (
        IOobject
        (
            "Ucyl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        Zero,
        dimVelocity,
        "calculated"
    );

    volVectorField URelcyl
    (
        IOobject
        (
            "URelcyl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        Zero,
        dimVelocity,
        "calculated"
    );

    vectorField radDirC = transform(planeDir, mesh_.C()) - transform(planeDir, origin_);
    radDirC /= mag(radDirC);
    vectorField tanDirC = axialDirection_ ^ radDirC;

    const vectorField& Uc = U.primitiveField();
    const vectorField& URelC = URel.primitiveField();
    vectorField& UcylC = Ucyl.primitiveFieldRef();
    vectorField& URelcylC = URelcyl.primitiveFieldRef();
    scalarField& UMerC = UMer.primitiveFieldRef();
    scalarField& URelMerC = URelMer.primitiveFieldRef();

    UcylC = ((axialDirection_&Uc)*vector(0.0,0.0,1.0)
           + (radDirC&Uc)*vector(1.0,0.0,0.0)
           + (tanDirC&Uc)*vector(0.0,1.0,0.0));

    URelcylC = ((axialDirection_&URelC)*vector(0.0,0.0,1.0)
              + (radDirC&URelC)*vector(1.0,0.0,0.0)
              + (tanDirC&URelC)*vector(0.0,1.0,0.0));

    UMerC = (mag((UcylC&vector(1.0,0.0,0.0))*vector(1.0,0.0,0.0)
          + (UcylC&vector(0.0,0.0,1.0))*vector(0,0,1)));

    URelMerC = (mag((URelcylC&vector(1.0,0.0,0.0))*vector(1.0,0.0,0.0)
             + (URelcylC&vector(0.0,0.0,1.0))*vector(0,0,1)));

    volScalarField h0("h0",h+0.5*magSqr(U));
    volScalarField p0("p0",p+0.5*rho*magSqr(U));
    volScalarField Rothalpy("Rothalpy", h + 0.5*magSqr(URel) - 0.5 * magSqr(Utr));

    Log << type() << " " << name() << " write:" << nl;

    if (!obr_.time().writeTime())
    {
        obr_.time().writeTimeDict();
    }

    for (const wordRe& objName : objectNames_)
    {
        Log << "    writing object " << objName << endl;
        if (objName == "Ucyl") Ucyl.write();
        else if (objName == "URelcyl") URelcyl.write();
        else if (objName == "URel") URel.write();
        else if (objName == "UMer") UMer.write();
        else if (objName == "URelMer") URelMer.write();
        else if (objName == "Utr") Utr.write();
        else if (objName == "h0") h0.write();
        else if (objName == "p0") p0.write();
        else if (objName == "Rothalpy") Rothalpy.write();
        else Info << "Invalid object: " << objName << endl;
    }

    return true;
}


// ************************************************************************* //
