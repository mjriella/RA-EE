    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", p.dimensions(), 0);
    dimensionedScalar pMax("pMax", p.dimensions(), 0);

    fieldBounds.lookup(p.name()) >> pMin.value() >> pMax.value();

    // Velocity bound
    dimensionedScalar U1Max("U1Max", U1.dimensions(), 0);
    dimensionedScalar U2Max("U2Max", U1.dimensions(), 0);

    fieldBounds.lookup(U1.name()) >> U1Max.value();
    fieldBounds.lookup(U2.name()) >> U2Max.value();

    dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
