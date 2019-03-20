if (turbulence1)
{
    if (mesh.changing())
    {
        y.correct();
    }

	// load turbulence constants
	#include "turbulenceConstants.H"

	// Drag coupling
	volScalarField dragCoup = 2*K/(alpha1*rho1);
    dragCoup.max(1e-8);

    // use N=6 so that f=0 at walls
    const dimensionedScalar N("N", dimless, 6.0);

    // original formulation - v2f
	volScalarField T = max(k1/epsilon1, 6*sqrt(nu2/epsilon2));
	volScalarField inter = pow3(nu2)/epsilon2;
	volScalarField L = CL*max(pow(k1, 1.5)/epsilon1, Ceta*pow(inter,0.25));

	volScalarField alpha
    (
        "v2f::alpha",
        1.0/T*((C1 - N)*v1 - 2.0/3.0*k1*(C1 - 1.0))
    );

    tmp<volScalarField> Ceps1 =
        1.4*(1.0 + 0.05*min(sqrt(k1/v1), scalar(100.0)));

    volScalarField G(2*nut1*magSqr(symm(fvc::grad(U1))));
    
    // Dissipation equation
    fvScalarMatrix epsEqn
    (
        fvm::ddt(epsilon1)
      + fvm::div(phi1, epsilon1)
      - fvm::laplacian
        (
            nut1/sigmaEps + nuKT, epsilon1,
            "laplacian(DepsilonpEff,epsilon1)"
        )
      ==
         Ceps1*G/T
       - fvm::Sp(Ceps2/T, epsilon1)
        + dragCoup*epsilonCoup
       - fvm::SuSp(dragCoup, epsilon1)
    );
    epsEqn.relax();
    epsEqn.solve();
	#include "../wallManipulation/epsilon1BC.H"
    #include "../wallManipulation/bounding/boundepsilon1.H"

    // Turbulent kinetic energy equation
    fvScalarMatrix k1Eqn
    (
        fvm::ddt(k1)
      + fvm::div(phi1, k1)
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

    // Relaxation function equation
    fvScalarMatrix f1Eqn
    (
      - fvm::laplacian(f1)
     ==
      - fvm::Sp(1.0/sqr(L), f1)
	  - 1.0/sqr(L)/k1*(alpha - C2*G)
    );
    f1Eqn.relax();
    f1Eqn.solve();
    #include "../wallManipulation/bounding/boundf1.H"

    // Turbulence stress normal to streamlines equation
    fvScalarMatrix v1Eqn
    (
        fvm::ddt(v1)
      + fvm::div(phi1, v1)
      - fvm::laplacian
        (
            nut1/sigmaK + nuKT, v1,
            "laplacian(DkEff,k1)"
        )
      ==
        min(k1*f1, -alpha + C2*G)
      - fvm::Sp(N*epsilon1/k1, v1)
      + dragCoup*vCoup
      - fvm::SuSp(dragCoup, v1)
    );
    v1Eqn.relax();
    v1Eqn.solve();
    #include "../wallManipulation/bounding/boundv1.H"

    //- Re-calculate turbulence viscosity
    nut1 = min(CmuKEps*sqr(k1)/epsilon1, Cmu*v1*T);
}

nuEff1 = nut1;
