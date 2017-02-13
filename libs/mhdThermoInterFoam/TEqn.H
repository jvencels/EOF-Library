{

    volScalarField kappaf = alpha1*rho1*cp1/Pr1*mixture.nuModel1().nu()
                              + alpha2*rho2*cp2/Pr2*mixture.nuModel2().nu();

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
