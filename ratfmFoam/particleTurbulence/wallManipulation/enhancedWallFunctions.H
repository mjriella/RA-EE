{
    labelList cellBoundaryFaceCount(epsilon1.size(), 0);

    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu.value(), 0.75);
    scalar kappa_ = kappa.value();
    scalar nu2_ = nu2.value();
    //scalar nu1_ = nu1.value();

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

                epsilon1[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    // Enhanced wall treatment has been added to save the 
    // turbulence at the wall

    //- Accumulate the wall face contributions to epsilon and G
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const scalarField& nut1w = nut1.boundaryField()[patchi];
            //const scalarField& nuKTw = nuKT.boundaryField()[patchi];

            scalarField magFaceGradU(mag(U1.boundaryField()[patchi].snGrad()));


            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                scalar yPlus = Cmu25*y[patchi][facei]*::sqrt(k2[faceCelli])/nu2_;

                scalar yPlusLam = 11.0;

                // For corner cells (with two boundary or more faces),
                // epsilon and G in the near-wall cell are calculated
                // as an average2

                cellBoundaryFaceCount[faceCelli]++;

                if (yPlus > yPlusLam)
                {
                    epsilon1[faceCelli] +=
                        Cmu75*::pow(k1[faceCelli], 1.5)
                        /(kappa_*y[patchi][facei]);

                    G[faceCelli] +=
                        (nut1w[facei] + nu2_)*magFaceGradU[facei]
                       *Cmu25*::sqrt(k1[faceCelli])
                        /(kappa_*y[patchi][facei]);
                }
                else
                {
                    // Note : Should check the viscosity used here
                    epsilon1[faceCelli] += 2.0*k1[faceCelli]*nu2_/
                                    sqr(y[patchi][facei]);
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

                epsilon1[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                G[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                cellBoundaryFaceCount[faceCelli] = 1;
            }
        }
    }

}
