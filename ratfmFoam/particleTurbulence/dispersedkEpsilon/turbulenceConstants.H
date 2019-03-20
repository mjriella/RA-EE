
    dimensionedScalar C1
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            k1Epsilon1Dict,
            1.44
        )
    );

    dimensionedScalar C2
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            k1Epsilon1Dict,
            1.9
        )
    );

    dimensionedScalar Cmu
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            k1Epsilon1Dict,
            0.09
        )
    );

    dimensionedScalar sigmaK
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            v2fDict,
            1.0
        )
    );

    dimensionedScalar sigmaEps
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            v2fDict,
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




