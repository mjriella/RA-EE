    IOdictionary turbulence1Properties
    (
        IOobject
        (
            "turbulence1Properties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary turbulenceModel1
    (
        turbulence1Properties.subDictPtr("turbulenceModel1")
    );

    Switch dispersedkEpsilon
    (
        turbulenceModel1.lookup("dispersedkEpsilon", false)
    );

    Switch dispersedv2f
    (
        turbulenceModel1.lookup("dispersedv2f", false)
    );

    Switch turbulence1
    (
        turbulence1Properties.lookup("turbulence")
    );

    // For coefficients of turbulence models
    dictionary k1Epsilon1Dict
    (
        turbulence1Properties.subDictPtr("k1Epsilon1Coeffs")
    );

    nearWallDist y(mesh);


    Info<< "Reading field k1\n" << endl;
    volScalarField k1
    (
        IOobject
        (
            "k1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon1\n" << endl;
    volScalarField epsilon1
    (
        IOobject
        (
            "epsilon1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field v2p\n" << endl;
    volScalarField v1
    (
        IOobject
        (
            "v1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField f1
    (
        IOobject
        (
            "f1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field nut1\n" << endl;
    volScalarField nut1
    (
        IOobject
        (
            "nut1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field nuEff1\n" << endl;
    volScalarField nuEff1
    (
        IOobject
        (
            "nuEff1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
    	dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    );




