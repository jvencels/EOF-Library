{
    // Effective thermal diffusivity: conductive + turbulent part
    volScalarField kappaf = kappa/cp + rho*turbulence->nut()/0.7;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(rhoPhi, T)
      - (
            fvm::laplacian(kappaf, T)
          + JH/cp
        )
    );

    TEqn.relax();
    TEqn.solve();
}
