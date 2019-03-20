{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar kappa_ = kappa.value();
	scalar E_ = E.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nut2.boundaryField()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // calculate yPlus
                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k2[faceCelli])/nu2.value();

                if (yPlus > 11.6)
                {
                    nutw[facei] = nu2.value()*(yPlus*kappa_/::log(E_*yPlus) - 1.0);
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}
