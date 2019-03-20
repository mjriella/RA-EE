    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar alphaSmall
    (
        "alphaSmall",
        dimless,
        1e-6
    );

    autoPtr<phaseModel> phase1 = phaseModel::New
    (
        mesh,
        transportProperties,
        "1"
    );

    volVectorField& U1 = phase1->U();
    surfaceScalarField& phi1 = phase1->phi();
    const dimensionedScalar& rho1 = phase1->rho();

    autoPtr<phaseModel> phase2 = phaseModel::New
    (
        mesh,
        transportProperties,
        "2"
    );

    volVectorField& U2 = phase2->U();
    surfaceScalarField& phi2 = phase2->phi();
    const dimensionedScalar& rho2 = phase2->rho();
    const dimensionedScalar& nu2 = phase2->nu();

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField alpha2
    (
        IOobject
        (
            "alpha2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(1) - alpha1
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1*U1 + alpha2*U2
    );

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(alpha1)*phi1 + fvc::interpolate(alpha2)*phi2
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        alpha1*rho1 + alpha2*rho2
    );

    IOdictionary interfacialProperties
    (
        IOobject
        (
            "interfacialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<dragModel> drag1 = dragModel::New
    (
        interfacialProperties,
        alpha1,
        phase1,
        phase2
    );

    kineticTheoryModel kineticTheory
    (
        phase1,
        U2,
        alpha1,
        drag1
    );

	// add lookup later
    volVectorField g_v
    (
        IOobject
        (
            "g_v",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet (0, 1, -2, 0, 0), vector(-9.8,0,0))
    );

    Info << "Creating field U1U2p\n" << endl;
    volVector7Field U1U2p
    (
        IOobject
        (
            "U1U2p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector7("zero", dimless, vector7::zero)
    );

    #include "createParticleTurbulence.H"
    #include "createFluidTurbulence.H"

    surfaceScalarField rAU1f
    (
        IOobject
        (
            "rAU1f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 0.0)
    );

    surfaceScalarField ppMagf
    (
        IOobject
        (
            "ppMagf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
	mesh.schemesDict().setFluxRequired(p.name());
