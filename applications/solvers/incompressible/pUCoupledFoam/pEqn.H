    surfaceScalarField presSource
    (
        "presSource",
        segregatedRhieChow
      ? pAlpha*fvc::flux(rAU*fvc::grad(p))
      : pAlpha*fvc::interpolate(rAU)
       *(fvc::interpolate(fvc::grad(p)) & mesh.Sf())
    );

    fvScalarMatrix pEqn
    (
        - pAlpha*fvm::laplacian(rAU, p)
     ==
        - fvc::div(presSource)
    );

    pEqn.setReference(pRefCell, pRefValue);

    UpEqn.insertEquation(pEqn, 3, 3);
