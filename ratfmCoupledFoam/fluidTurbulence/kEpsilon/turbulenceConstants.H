
    dimensionedScalar C1
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            kEpsilonDict,
            1.44
        )
    );

    dimensionedScalar C2
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            kEpsilonDict,
            1.92
        )
    );

    dimensionedScalar Cmu
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            kEpsilonDict,
            0.09
        )
    );

    dimensionedScalar sigmaK
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            kEpsilonDict,
            1.0
        )
    );

    dimensionedScalar sigmaEps
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            kEpsilonDict,
            1.3
        )
    );

    dictionary wallFunctionDict
    (
        turbulence2Properties.subDictPtr("wallFunctionCoeffs")
    );

    dimensionedScalar kappa
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            wallFunctionDict,
            0.41
        )
    );

    dimensionedScalar E
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "E",
            wallFunctionDict,
            9.8
        )
    );




