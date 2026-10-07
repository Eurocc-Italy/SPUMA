/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2025-2026 Cineca
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

Application
    pUCoupledFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state incompressible solver with implicit pressure-velocity
    coupling in a vector4 block system.

    The Rhie-Chow flux is selected with the fvSolution entry rhieChow:
    coupled (default, interpolate(rAU)*interpolate(grad(p))) or segregated
    (flux(rAU*grad(p)) with rAU of the relaxed UEqn, as simpleFoam).

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvBlockMatrix.H"
#include "fvBlockMatrixOperators.H"
#include "vector4.H"
#include "tensor4.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state incompressible pressure-velocity block solver."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "initDevice.H"
    #include "createMemoryPool.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    const scalar pAlpha = mesh.solutionDict().getOrDefault("pAlpha", 1.0);

    const word rhieChow
    (
        mesh.solutionDict().getOrDefault<word>("rhieChow", "coupled")
    );

    if (rhieChow != "coupled" && rhieChow != "segregated")
    {
        FatalIOErrorInFunction(mesh.solutionDict())
            << "Unknown rhieChow " << rhieChow
            << ", valid options: coupled segregated" << exit(FatalIOError);
    }

    const bool segregatedRhieChow = (rhieChow == "segregated");

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvBlockMatrix<vector4, tensor4> UpEqn(mesh);

        #include "UEqn.H"
        #include "pEqn.H"
        #include "couplingTerms.H"

        UpEqn.solve(U, 0, p, 3);

        phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

        #include "continuityErrs.H"

        fvOptions.correct(U);
        fvOptions.correct(p);

        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

        #include "poolOccupancy.H"
        runTime.printExecutionTime(Info);
    }

    #include "poolMaxOccupancy.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
