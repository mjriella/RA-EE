{
    labelList cellBoundaryFaceCount(epsilon2.size(), 0);

    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu.value(), 0.75);
    scalar kappa_ = kappa.value();
    scalar nu2_ = nu2.value();

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                epsilon2[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to epsilon and G
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const scalarField& nut2w = nut2.boundaryField()[patchi];

            scalarField magFaceGradU(mag(U2.boundaryField()[patchi].snGrad()));


            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]
                    *::sqrt(k2[faceCelli])
                    /nu2.value();

                cellBoundaryFaceCount[faceCelli]++;

                epsilon2[faceCelli] +=
                     Cmu75*::pow(k2[faceCelli], 1.5)
                    /(kappa_*y[patchi][facei]);

                if (yPlus > 11.6)
                {
                G[faceCelli] +=
                    sqr((nut2w[facei]+nu2_)*magFaceGradU[facei])/
                   (Cmu25*::sqrt(k2[faceCelli])*kappa_*y[patchi][facei]);
                }
				else
				{
				G[faceCelli] = 0.0;
				}	

            }
        }
    }


    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon2[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                G[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                cellBoundaryFaceCount[faceCelli] = 1;
            }
        }
    }

}
