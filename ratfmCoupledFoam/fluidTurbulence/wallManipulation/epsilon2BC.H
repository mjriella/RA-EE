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
    		dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), turbulence2Properties.lookup("nu"))
		)
);

volScalarField& nu2 = nu_temp();

{
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];
		const scalarField& nu2f = nu2.boundaryField()[patchi];

        if (isA<wallFvPatch>(p))
        {
			forAll(p, facei)
			{
				label faceCelli = p.faceCells()[facei];
				
				epsilon2[faceCelli] =
					2*k2[faceCelli]*nu2f[facei]/
							sqr(y[patchi][facei]);
			}
        }
    }
}


