{
    volScalarField kappaf = alpha1f*kappa1/cp1 + (scalar(1)-alpha1f)*kappa2/cp2 + rho*turbulence->nut()/0.7;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - (
            fvm::laplacian(kappaf, T)
          + JH*(alpha1f/cp1 + (scalar(1)-alpha1f)/cp2)
        )
    );

    TEqn.relax();
    TEqn.solve();
}
