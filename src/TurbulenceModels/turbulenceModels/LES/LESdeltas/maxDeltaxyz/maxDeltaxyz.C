/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "maxDeltaxyz.H"
#include "executors.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(maxDeltaxyz, 0);
    addToRunTimeSelectionTable(LESdelta, maxDeltaxyz, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyz::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    label nD = mesh.nGeometricD();

    const auto cellFaceStart = mesh.cellsFaceStart().cbegin();
    const auto cellFaces = mesh.cellsFaces().cbegin();
    const vectorField& cellC = mesh.cellCentres();
    const vectorField& faceC = mesh.faceCentres();
    const vectorField faceN(mesh.faceAreas()/mag(mesh.faceAreas()));
    scalarField hmax(mesh.nCells());

    const auto cellCp = cellC.cbegin();
    const auto faceCp = faceC.cbegin();
    const auto faceNp = faceN.cbegin();
    auto hmaxp = hmax.begin();
    const scalar deltaCoeff = deltaCoeff_;
    foamExecutor exec;
    auto calculate = [=](label celli)
    {
        scalar deltaMaxTmp = 0.0;
        const point& cc = cellCp[celli];

        for (label i = cellFaceStart[celli]; i < cellFaceStart[celli+1]; ++i)
        {
            label facei = cellFaces[i];
            const point& fc = faceCp[facei];
            const vector& n = faceNp[facei];

            scalar tmp = mag(n & (fc - cc));
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        hmaxp[celli] = deltaCoeff*deltaMaxTmp;
    };
    exec.parallelFor(calculate, mesh.nCells());

    if (nD == 3)
    {
        delta_.primitiveFieldRef() = hmax;
    }
    else if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;

        delta_.primitiveFieldRef() = hmax;
    }
    else
    {
        FatalErrorInFunction
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::maxDeltaxyz::maxDeltaxyz
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            2
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::maxDeltaxyz::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);

    calcDelta();
}


void Foam::LESModels::maxDeltaxyz::correct()
{
    if (turbulenceModel_.mesh().changing())
    {
        calcDelta();
    }
}


// ************************************************************************* //
