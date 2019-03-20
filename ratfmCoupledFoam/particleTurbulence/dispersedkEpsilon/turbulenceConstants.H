
    dimensionedScalar C1
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            dispersedkEpsilonDict,
            1.44
        )
    );

    dimensionedScalar C2
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            dispersedkEpsilonDict,
            1.9
        )
    );

    dimensionedScalar Cmu
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            dispersedkEpsilonDict,
            0.09
        )
    );

    dimensionedScalar sigmaK
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            dispersedkEpsilonDict,
            1.0
        )
    );

    dimensionedScalar sigmaEps
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            dispersedkEpsilonDict,
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




