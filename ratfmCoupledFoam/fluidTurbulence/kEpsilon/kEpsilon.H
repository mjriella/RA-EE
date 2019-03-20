if (turbulence2)
{
    if (mesh.changing())
    {
        y.correct();
    }

	// Add constants 
	#include "turbulenceConstants.H"

	// Drag coupling - divide by phase 
	volScalarField dragCoup = 2*K/((1-alpha1)*rho2);

    dragCoup.max(1e-8);

    // Construct production term
    volScalarField G(2*nut2*magSqr(symm(fvc::grad(U2))));

    #include "../wallManipulation/enhancedWallFunctions.H"

    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon2)
      + fvm::div(phi2, epsilon2)
      - fvm::Sp(fvc::div(phi2), epsilon2)
      - fvm::laplacian(nut2/sigmaEps + nu2, epsilon2)
     ==
        C1*G*epsilon2/k2
      - fvm::Sp(C2*epsilon2/k2, epsilon2)
      + dragCoup*epsilonCoup 
      - fvm::SuSp(dragCoup, epsilon2)
    );

    #include "../wallManipulation/wallDissipation.H"

    epsEqn.relax();
    epsEqn.solve();

    #include "../wallManipulation/bounding/boundepsilon2.H"


    // Turbulent k2inetic energy equation
    fvScalarMatrix k2Eqn
    (
        fvm::ddt(k2)
      + fvm::div(phi2, k2)
      - fvm::Sp(fvc::div(phi2), k2)
      - fvm::laplacian(nut2/sigmaK + nu2, k2)
     ==
        G
      - fvm::Sp(epsilon2/k2, k2)
      + dragCoup*kCoup
      - fvm::SuSp(dragCoup, k2)
    );
    k2Eqn.relax();
    k2Eqn.solve();

    #include "../wallManipulation/bounding/boundk2.H"
    //- Re-calculate turbulence viscosity
    nut2 = Cmu*sqr(k2)/epsilon2; 
    #include "../wallManipulation/wallViscosity.H"
}
nuEff2 = nut2 + nu2;
