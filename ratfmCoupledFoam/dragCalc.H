volVectorField Ur(U1 - U2);
volScalarField magUr(mag(Ur));

volScalarField Ka(drag1->K(magUr));
volScalarField K("K", Ka);

// Remove lift and drag at fixed-flux boundaries
forAll(phi1.boundaryField(), patchi)
{
    if (isA<fixedValueFvsPatchScalarField>(phi1.boundaryField()[patchi]))
    {
        K.boundaryField()[patchi] = 0.0;
    }
}

