{
    surfaceScalarField alpha1f(fvc::interpolate(alpha1));
    surfaceScalarField alpha2f(scalar(1) - alpha1f);

    volScalarField rAU1(1.0/U1Eqn.A());
    volScalarField rAU2(1.0/U2Eqn.A());

    rAU1f= 1.0/fvc::interpolate(U1Eqn.A());
    surfaceScalarField rAU2f(1.0/fvc::interpolate(U2Eqn.A()));

    volVectorField HbyA1("HbyA1", U1);
    HbyA1 = rAU1*U1Eqn.H();

    volVectorField HbyA2("HbyA2", U2);
    HbyA2 = rAU2*U2Eqn.H();

    surfaceScalarField ppDrag("ppDrag", 0.0*phi1);

    if (g0.value() > 0.0)
    {
        ppDrag -= ppMagf*fvc::snGrad(alpha1)*mesh.magSf();
    }

    if (kineticTheory.on())
    {
        ppDrag -= rAU1f*fvc::snGrad(kineticTheory.pa()/rho1)*mesh.magSf();
    }

    // Turbulent diffusion
    volScalarField Scfp = sqrt(k2/k1);
    volScalarField D = nut2/(alpha1*alpha2*Scfp);
    
    surfaceScalarField Df1(rAU1f*fvc::interpolate(D*drag12)); 
    surfaceScalarField Df2(rAU2f*fvc::interpolate(D*drag21));

    // Phase-fraction gradient (face) - transfered from turbDisp model
    surfaceScalarField snGradAlpha(fvc::snGrad(alpha1)*mesh.magSf());

    surfaceScalarField phiD1 = -Df1*snGradAlpha;
    surfaceScalarField phiD2 =  Df2*snGradAlpha;

    // Phase-1 predicted flux
    surfaceScalarField phiHbyA1
    (
        "phiHbyA1",
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, U1, phi1) 
    );

    // Phase-2 predicted flux
    surfaceScalarField phiHbyA2
    (
        "phiHbyA2",
        (fvc::interpolate(HbyA2) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU2, U2, phi2)
    );

    phi = alpha1f*phiHbyA1 + alpha2f*phiHbyA2;

    phiHbyA1 +=
    (
        fvc::interpolate(drag21*rAU1)*phi2
      + ppDrag
      + rAU1f*(g & mesh.Sf())
	  + phiD1
    );

    phiHbyA2 +=
    (
        fvc::interpolate(drag21*rAU2)*phi1
      + rAU2f*(g & mesh.Sf())
      + phiD2
    );
   
    // Mean predicted flux for continuity
    surfaceScalarField phiHbyA("phiHbyA", alpha1f*phiHbyA1 + alpha2f*phiHbyA2);

    HbyA1 += drag12*rAU1*U2;
    HbyA2 += drag21*rAU2*U1;

    // Pressure diffusivity 
    surfaceScalarField Dp
    (
        "Dp",
        (
            alpha1f*rAU1f/rho1 
          + alpha2f*rAU2f/rho2
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        // solve p with predicted fluxes
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        maxResidualp = pEqn.solve
        (
            mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
        ).initialResidual();


        // Correct fluxes and velocities on last iteration
        if (pimple.finalNonOrthogonalIter())
        {
            surfaceScalarField SfGradp(pEqn.flux()/Dp);

            // correct fluxes with mass conservative field
            phi1 = phiHbyA1 - rAU1f*SfGradp/rho1;
            phi2 = phiHbyA2 - rAU2f*SfGradp/rho2;
            phi = alpha1f*phi1 + alpha2f*phi2;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

	    	// correct velocities with mass conservative field
            U1 = HbyA1
               + fvc::reconstruct
                 (
                     ppDrag
                   + rAU1f*((g & mesh.Sf()) - SfGradp/rho1)
                   + phiD1
                 );
            U1.correctBoundaryConditions();

            U2 = HbyA2
               + fvc::reconstruct
                 (
                     rAU2f*((g & mesh.Sf()) - SfGradp/rho2)
                   + phiD2
                 );
            U2.correctBoundaryConditions();

            U = alpha1*U1 + alpha2*U2;
        }
    }
}

#include "continuityErrs.H"
