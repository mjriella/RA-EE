surfaceScalarField alpha1f(fvc::interpolate(alpha1));
surfaceScalarField alpha2f(scalar(1) - alpha1f);

volScalarField rAU1(1.0/U1Eqn.A());
volScalarField rAU2(1.0/U2Eqn.A());

rAU1f= 1.0/fvc::interpolate(U1Eqn.A());
surfaceScalarField rAU2f(1.0/fvc::interpolate(U2Eqn.A()));

// Define coefficient ratios in U1
volScalarField diagU1n = U1Eqn.A()/U1Eqn.A();
// drag in 1
volScalarField thetaD1 = drag12*rAU1;
surfaceScalarField thetaD1f = fvc::interpolate(thetaD1);

// Define coefficient ratios in U2
volScalarField diagU2n = U2Eqn.A()/U2Eqn.A();
// drag in 2
volScalarField thetaD2 = drag21*rAU2;
surfaceScalarField thetaD2f = fvc::interpolate(thetaD2);

// Time in U1 
fvVectorMatrix U1T = fvm::ddt(U1);
volScalarField thetaT1 = U1T.A()*rAU1;
surfaceScalarField thetaT1f = fvc::interpolate(thetaT1);
// Time in U2 
fvVectorMatrix U2T = fvm::ddt(U2);
volScalarField thetaT2 = U2T.A()*rAU2;
surfaceScalarField thetaT2f = fvc::interpolate(thetaT2);


// division ratios
volScalarField oneOverRatio1 = diagU1n + thetaT1 + thetaD1;
surfaceScalarField oneOverRatio1f = fvc::interpolate(oneOverRatio1);
volScalarField oneOverRatio2 = diagU2n + thetaT2 + thetaD2;
surfaceScalarField oneOverRatio2f = fvc::interpolate(oneOverRatio2);

// explicit drag in pEqn
surfaceScalarField phiDrag1 = (thetaD1f*-phi2)/oneOverRatio1f;
surfaceScalarField phiDrag2 = (thetaD2f*-phi1)/oneOverRatio2f;

surfaceScalarField phiDragt1("phiDragt1", alpha1f*phiDrag1);
surfaceScalarField phiDragt2("phiDragt2", alpha2f*phiDrag2);

// corrective drag - implicit div
// discretise opposite contribution 
surfaceScalarField correctionU1 = alpha1f 
                            - ((alpha2f*thetaD2f)/oneOverRatio2f);
surfaceScalarField correctionU2 = alpha2f 
                            - ((alpha1f*thetaD1f)/oneOverRatio1f);

// exp time 
surfaceScalarField phit1 = (thetaT1f*((fvc::interpolate(U1.oldTime()) & mesh.Sf()) 
                                    - phi1.oldTime()))/oneOverRatio1f;
surfaceScalarField phit2 = (thetaT2f*((fvc::interpolate(U2.oldTime()) & mesh.Sf()) 
                                    -phi2.oldTime()))/oneOverRatio2f;

surfaceScalarField phiT1("phit1", alpha1f*phit1);
surfaceScalarField phiT2("phit2", alpha2f*phit2);

// pressure diffusivity
volScalarField Dp
(
   "Dp",
    (
       	alpha1/(rho1*oneOverRatio1)*rAU1
      + alpha2/(rho2*oneOverRatio2)*rAU2
    )
);

surfaceScalarField Dp_f("Dp_f", fvc::interpolate(Dp));

// explicit interpolated flux in Rhie-Chow - 1/s 
surfaceScalarField pressureSource
(
    "pressureSource",
    Dp_f*(fvc::interpolate(fvc::grad(p)) & mesh.Sf())
);

// correct flux
surfaceScalarField SourceSfp(pressureSource/Dp_f);

// pEqn is here now
fvScalarMatrix pEqn
(	
    - fvm::laplacian(Dp, p)
    ==
    - fvc::div(pressureSource)
    + fvc::div(phiDragt1)
    + fvc::div(phiDragt2)
    + fvc::div(phiT1)
    + fvc::div(phiT2)
);

U1U2pEqn.insertEquation(6, pEqn);
