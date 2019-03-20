{
    // Bound the pressure
    dimensionedScalar p1 = min(p);
    dimensionedScalar p2 = max(p);

    if (p1 < pMin || p2 > pMax)
    {
        Info<< "p: " << p1.value() << " " << p2.value()
            << ".  Bounding." << endl;

        p.max(pMin);
        p.min(pMax);
        p.correctBoundaryConditions();
    }

    // Bound U1
    volScalarField magU1 = mag(U1);
    dimensionedScalar U1Lim = max(magU1);

    if (U1Lim > U1Max)
    {
        Info<< "U1: " << U1Lim.value() << ".  Bounding." << endl;

        volScalarField U1limiter = pos(magU1 - U1Max)*U1Max/(magU1 + smallU)
            + neg(magU1 - U1Max);
        U1limiter.max(scalar(0));
        U1limiter.min(scalar(1));

        U1 *= U1limiter;
        U1.correctBoundaryConditions();
    }

    // Bound U2
    volScalarField magU2 = mag(U2);
    dimensionedScalar U2Lim = max(magU2);

    if (U2Lim > U2Max)
    {
        Info<< "U2: " << U2Lim.value() << ".  Bounding." << endl;

        volScalarField U2limiter = pos(magU2 - U2Max)*U2Max/(magU2 + smallU)
            + neg(magU2 - U2Max);
        U2limiter.max(scalar(0));
        U2limiter.min(scalar(1));

        U2 *= U2limiter;
        U2.correctBoundaryConditions();
    }
}
