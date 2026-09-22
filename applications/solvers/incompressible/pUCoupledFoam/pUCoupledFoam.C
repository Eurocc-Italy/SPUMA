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
    Steady-state incompressible solver scaffold for a directly assembled
    vector4 pressure-velocity block system.

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

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Work in progress: assemble the vector4/tensor4 coupled system
        // to exercise the new insertion functions before adding the solve.
        fvBlockMatrix<vector4, tensor4> UpEqn(mesh);

        // Assemble the momentum equation
        fvVectorMatrix UEqn
        (
                fvm::div(phi, U)
            + turbulence->divDevReff(U)
            - fvOptions(U)
        );

        // Store the reciprocal momentum diagonal before relaxation
        volScalarField rAU(1.0/UEqn.A());

        // Relax and insert the momentum block
        UEqn.relax();
        UpEqn.insertMomentumEquation(UEqn);

        // Build the Rhie-Chow pressure-correction source
        surfaceScalarField presSource
        (
                "presSource",
                pAlpha*fvc::interpolate(rAU)
             *(fvc::interpolate(fvc::grad(p)) & mesh.Sf())
        );

        // Assemble the pressure equation
        fvScalarMatrix pEqn
        (
            - pAlpha*fvm::laplacian(rAU, p)
         ==
            - fvc::div(presSource)
        );

        // Set the pressure reference and insert the pressure block
        pEqn.setReference(pRefCell, pRefValue);
        UpEqn.insertPressureEquation(pEqn);

        // Assemble the off-diagonal pressure-velocity coupling blocks
        fvBlockMatrix<vector, vector> dUByp(fvm::grad(p));
        fvBlockMatrix<scalar, vector> dpByU(fvm::div(U));

        UpEqn.insertVelocityPressureCoupling(dUByp);
        UpEqn.insertPressureVelocityCoupling(dpByU);

        // Assembly-only test: U and p are intentionally not updated yet.
    }

    #include "poolMaxOccupancy.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
