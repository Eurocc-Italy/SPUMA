    fvBlockMatrix<vector, vector> dUByp(fvm::grad(p));
    fvBlockMatrix<scalar, vector> dpByU(fvm::div(U));

    UpEqn.insertEquation(dUByp, 0, 3);
    UpEqn.insertEquation(dpByU, 3, 0);
