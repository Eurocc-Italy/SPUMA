    fvVectorMatrix UEqn
    (
            fvm::div(phi, U)
        + turbulence->divDevReff(U)
        - fvOptions(U)
    );

    volScalarField rAU(1.0/UEqn.A());

    if (segregatedRhieChow)
    {
        // Relaxation only enters rAU; the block system is not relaxed
        fvVectorMatrix UEqnRelaxed(UEqn);
        UEqnRelaxed.relax();
        rAU = 1.0/UEqnRelaxed.A();
    }
    else
    {
        UEqn.relax();
    }

    UpEqn.insertEquation(UEqn, 0, 0);
