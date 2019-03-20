if (turbulence2)
{
    if (mesh.changing())
    {
        y.correct();
    }

	// Load constants 
	#include "turbulenceConstants.H"

	volScalarField dragCoup = 2*K/((1-alpha1)*rho2);
    dragCoup.max(1e-8);

	// Coupling in vel-pres model
	volScalarField vCoup = kCoup;
	if(dispersedkEpsilon)
	{
		vCoup = sqrt((2/3)*v2*k1);
	}	
	else
	{
		vCoup = sqrt(v2*v1);
	}

    // use N=6 so that f=0 at walls
    const dimensionedScalar N("N", dimless, 6.0);

	volScalarField T = max(k2/epsilon2, 6*sqrt(nu2/epsilon2));
	volScalarField inter = pow3(nu2)/epsilon2;
	volScalarField L = CL*max(pow(k2, 1.5)/epsilon2, Ceta*pow(inter, 0.25));
	volScalarField alpha
    (
        "v2f::alpha",
        1.0/T*((C1 - N)*v2 - 2.0/3.0*k2*(C1 - 1.0))
    );

    tmp<volScalarField> Ceps1 =
        1.4*(1.0 + 0.05*min(sqrt(k2/v2), scalar(100.0)));

    // Construct production term
    volScalarField G(2*nut2*magSqr(symm(fvc::grad(U2))));

    // Wall BC for epsilon
	#include "../wallManipulation/epsilon2BC.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon2)
      + fvm::div(phi2, epsilon2)
	  - fvm::Sp(fvc::div(phi2), epsilon2)
      - fvm::laplacian(nut2/sigmaEps + nu2, epsilon2)
     ==
        Ceps1*G/T // 1.6 
      - fvm::Sp(Ceps2/T, epsilon2) 
      + dragCoup*epsilonCoup
      - fvm::SuSp(dragCoup, epsilon2)
    );
    epsEqn().relax();
    epsEqn().solve();
    #include "../wallManipulation/bounding/boundepsilon2.H"

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> k2Eqn
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
    k2Eqn().relax();
    k2Eqn().solve();
    #include "../wallManipulation/bounding/boundk2.H"

    // Relaxation function equation
    fvScalarMatrix f2Eqn
    (
      - fvm::laplacian(f2)
     ==
      - fvm::Sp(1.0/sqr(L), f2)
	  - 1.0/sqr(L)/k2*(alpha - C2*G)
    );
    f2Eqn.relax();
    f2Eqn.solve();
    #include "../wallManipulation/bounding/boundf.H"

    // Turbulence stress normal to streamlines equation
    fvScalarMatrix v2Eqn
    (
        fvm::ddt(v2)
      + fvm::div(phi2, v2)
      - fvm::laplacian(nut2/sigmaK + nu2, v2)
      ==
        min(k2*f2, -alpha + C2*G)
      - fvm::Sp(N*epsilon2/k2, v2)
      + dragCoup*vCoup
      - fvm::SuSp(dragCoup, v2)
    );
    v2Eqn.relax();
    v2Eqn.solve();
    
    #include "../wallManipulation/bounding/boundv2.H"

    //- Re-calculate turbulence viscosity
    nut2 = min(CmuKEps*sqr(k2)/epsilon2, Cmu*v2*T); // Davidson correct
}

nuEff2 = nut2 + nu2;
