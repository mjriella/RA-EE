    IOdictionary turbulence2Properties
    (
        IOobject
        (
            "turbulence2Properties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dictionary turbulenceModel2
    (
        turbulence2Properties.subDictPtr("turbulenceModel2")
    );

    Switch turbulence2
    (
        turbulence2Properties.lookup("turbulence")
    );

    Switch kEpsilon
    (
        turbulenceModel2.lookup("kEpsilon", false)
    );

    Switch v2f
    (
        turbulenceModel2.lookup("v2f", false)
    );

    dictionary v2fDict
    (
        turbulence2Properties.subDictPtr("v2fDict")
    );

    Info<< "Reading field k2\n" << endl;
    volScalarField k2
    (
        IOobject
        (
            "k2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon2\n" << endl;
    volScalarField epsilon2
    (
        IOobject
        (
            "epsilon2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field v2\n" << endl;
    volScalarField v2
    (
        IOobject
        (
            "v2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field f2\n" << endl;
    volScalarField f2
    (
        IOobject
        (
            "f2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field nut2\n" << endl;
    volScalarField nut2
    (
        IOobject
        (
            "nut2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
		mesh
    );

    Info<< "Calculating field nuEff2\n" << endl;
    volScalarField nuEff2
    (
        IOobject
        (
            "nuEff2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
    	dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    );

    volScalarField kCoup
    (
        IOobject
        (
            "kCoup",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
    	dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
    );

    volScalarField vCoup
    (
        IOobject
        (
            "vCoup",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
    	dimensionedScalar("zero", dimensionSet(0, 2, -2, 0, 0), 0.0)
    );

    volScalarField epsilonCoup
    (
        IOobject
        (
            "epsilonCoup",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
		mesh,
    	dimensionedScalar("zero", dimensionSet(0, 2, -3, 0, 0), 0.0)
    );

