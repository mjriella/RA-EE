volScalarField rho1Field
(
    IOobject
    (
        "rho1Field",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rho1
);

volScalarField Scfp = sqrt(k/kp);
volScalarField alpha2 = 1-alpha1;

volScalarField beta = alpha1.mesh().lookupObject<volScalarField> ("beta");
// particle relaxation time (s)
volScalarField tau_d = rho1Field/
                        max(beta, betaSmall);

volVectorField dispersion("dispersion", nut2*fvc::grad(alpha1)/
            (
             max(alpha1, scalar(1e-6))*
             max(alpha2, scalar(1e-6))*
             Scfp
            ));
//Info << "disppppp" << max(dispersion) << endl;
//Info << "tau_d" << max(tau_d) << endl;


volVectorField Urel = U1 -U2;
volScalarField prodSlip = -dispersion & Urel/tau_d;

volVectorField buoyancy = nut2*fvc::grad(alpha1)/Scfp;

volScalarField prodBuoy = -buoyancy & g;

Info << "Slip = " << max(prodSlip).value() <<" min = " 
            << min(prodSlip).value() << endl; 
Info << "Bouy = " << max(prodBuoy).value() << " min = "
            <<min(prodBuoy).value() << endl; 
