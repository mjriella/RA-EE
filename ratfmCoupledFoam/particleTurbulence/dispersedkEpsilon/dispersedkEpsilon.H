if (turbulence1)
{
    if (mesh.changing())
    {
        y.correct();
    }
	
	#include "turbulenceConstants.H"

	// Drag coupling
	volScalarField dragCoup = 2*K/(alpha1*rho1);
    dragCoup.max(1e-8);
    volScalarField G(2*nut1*magSqr(symm(fvc::grad(U1))));

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon1)
      + fvm::div(phi1, epsilon1)
      + fvm::Sp(-fvc::div(phi1), epsilon1)
      - fvm::laplacian
        (
            nut1/sigmaEps + nuKT, epsilon1,
            "laplacian(Depsilon1Eff,epsilon1)"
        )
      ==
         C1*G*epsilon1/k1
      - fvm::Sp(C2*epsilon1/k1, epsilon1)
      + dragCoup*epsilonCoup
      - fvm::SuSp(dragCoup, epsilon1)
    );

    epsEqn.relax();
    epsEqn.solve();

    #include "../wallManipulation/bounding/boundepsilon1.H"

    // Turbulent kinetic energy equation
    fvScalarMatrix k1Eqn
    (
        fvm::ddt(k1)
      + fvm::div(phi1, k1)
      + fvm::Sp(-fvc::div(phi1), k1)
      - fvm::laplacian
        (
            nut1/sigmaK + nuKT, k1,
            "laplacian(DkEff,k1)"
        )
      ==
        G
      - fvm::Sp(epsilon1/k1, k1)
      + dragCoup*kCoup
      - fvm::SuSp(dragCoup, k1)
    );

    k1Eqn.relax();
    k1Eqn.solve();

    #include "../wallManipulation/bounding/boundk1.H"

    //- Re-calculate turbulence viscosity
    nut1 = Cmu*sqr(k1)/epsilon1;
    nut1.correctBoundaryConditions();
}
nuEff1 = nut1;
