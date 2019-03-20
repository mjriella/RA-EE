// Set epsilon wall value to 2*k*nu/y^2
// ref is durbin 2011

tmp<volScalarField> nu_temp
	(
		new volScalarField
		(
    		IOobject
    		(
        		"nu_temp",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
    		),
    		mesh,
			nu2
		)
);

volScalarField& nu2 = nu_temp();

{

    labelList cellBoundaryFaceCount(epsilon1.size(), 0);

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall epsilon fields to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon1[faceCelli] = 0.0;
            }
        }
    }

    //const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];
		const scalarField& nu2f = nu2.boundaryField()[patchi];

        if (isA<wallFvPatch>(p))
        {
			forAll(p, facei)
			{
				label faceCelli = p.faceCells()[facei];
				
				epsilon1[faceCelli] =
					2*k1[faceCelli]*nu2f[facei]/
							sqr(y[patchi][facei]);
			}
        }
    }
}


