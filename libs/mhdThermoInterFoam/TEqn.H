{
    volScalarField kappaf = alpha1f*kappa1 + (scalar(1)-alpha1f)*kappa2;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoPhiCp, T)
      - (
            fvm::laplacian(kappaf, T)
          + JH
        )
    );

    TEqn.relax();
    TEqn.solve();
}
