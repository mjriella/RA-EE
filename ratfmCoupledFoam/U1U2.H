fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, U2.dimensions()*dimVol/dimTime);

// create drag in both phases.
volScalarField drag12 = (alpha2*K)/rho1;
volScalarField drag21 = (alpha1*K)/rho2;

// Turbulent diffusion.
volScalarField Scfp = sqrt(k2/k1);
volScalarField D = nut2/(alpha1*alpha2*Scfp);
    
volScalarField Df1(D*drag12); 
volScalarField Df2(D*drag21);

// Phase-fraction gradient. 
volVectorField gradAlpha(fvc::grad(alpha1));

volVectorField D1 = -Df1*gradAlpha;
volVectorField D2 =  Df2*gradAlpha;

// Particle pressure.
volVectorField vPp = fvc::grad(kineticTheory.pa())/rho1;

// Discretise U1.
volTensorField gradU1T = fvc::grad(U1)().T();

if(!turbulence1)
{
    nuEff1 = nu2;
}

if (kineticTheory.on())
{
    kineticTheory.solve();
    nuEff1 += kineticTheory.mu1()/rho1;
}

volTensorField Rc1
(
    "Rc1",
    (((2.0/3.0)*I)*nuEff1)*tr(gradU1T) - nuEff1*gradU1T
);

surfaceScalarField phiR1
(
    -fvc::interpolate(nuEff1)*mesh.magSf()*fvc::snGrad(alpha1)
    /fvc::interpolate(alpha1 + alphaSmall)
);

U1Eqn =
(
    fvm::ddt(U1)
  + fvm::div(phi1, U1, "div(phi1,U1)")
  - fvm::Sp(fvc::div(phi1), U1)
  - fvm::laplacian(nuEff1, U1) 
  + fvc::div(Rc1)
  + fvm::div(phiR1, U1, "div(phi1,U1)") 
  - fvm::Sp(fvc::div(phiR1), U1)
  + (fvc::grad(alpha1)/(fvc::average(alpha1) + alphaSmall) & Rc1)
  + fvm::Sp(drag12, U1)
 ==
   drag12*U2 
  + D1
  + g
  + vPp
);

U1Eqn.relax();
U1U2pEqn.insertEquation(0, U1Eqn);

// Discretise U2.
if(!turbulence2)
{
    nuEff2 = nu2;
}

volTensorField gradU2T = fvc::grad(U2)().T();
volTensorField Rc2
(
    "Rc2",
    (((2.0/3.0)*I)*nuEff2)*tr(gradU2T) - nuEff2*gradU2T
);

surfaceScalarField phiR2
(
    -fvc::interpolate(nuEff2)*mesh.magSf()*fvc::snGrad(alpha2)
    /fvc::interpolate(alpha2 + alphaSmall)
);

U2Eqn =
(
    fvm::ddt(U2)
  + fvm::div(phi2, U2, "div(phi2,U2)")
  - fvm::Sp(fvc::div(phi2), U2)
  - fvm::laplacian(nuEff2, U2) 
  + fvc::div(Rc2)
  + fvm::div(phiR2, U2, "div(phi2,U2)") 
  - fvm::Sp(fvc::div(phiR2), U2)
  + (fvc::grad(alpha2)/(fvc::average(alpha2) + alphaSmall) & Rc2)
  + fvm::Sp(drag21, U2)
 ==
   drag21*U1
  + D2
  + g
);
U2Eqn.relax();
U1U2pEqn.insertEquation(3, U2Eqn);
    
