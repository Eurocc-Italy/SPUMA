    fvVectorMatrix UEqn
    (
            fvm::div(phi, U)
        + turbulence->divDevReff(U)
        - fvOptions(U)
    );

    // Reciprocal momentum diagonal before relaxation
    volScalarField rAU(1.0/UEqn.A());

    UEqn.relax();

    UpEqn.insertEquation(UEqn, 0, 0);
