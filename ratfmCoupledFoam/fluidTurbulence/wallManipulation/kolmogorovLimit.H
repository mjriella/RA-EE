// convert dim scalar nu to volScalarField
// Manipulate boundaries
tmp<volScalarField> nu_field
	(
		new volScalarField
		(
    		IOobject
    		(
        		"nu_field",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
    		),
    		mesh,
    		dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), turbulence2Properties.lookup("nu"))
		)
);

volScalarField& nu_lam = nu_field();



// wall unit as defined by nu2/sqrt(taufr/rho2)
volScalarField ystar
(
    IOobject
    (
       "ystar",
        mesh.time().constant(),
        mesh
    ),
    mesh,
    dimensionedScalar("ystar", dimLength, 1.0)
);

const fvPatchList& patches = mesh.boundary();
forAll(patches, patchi)
{
    if (isA<wallFvPatch>(patches[patchi]))
    {
        const fvPatchVectorField& Uw = U2.boundaryField()[patchi];
        const scalarField& nuw = nu_lam.boundaryField()[patchi];

        ystar.boundaryField()[patchi] =
            nuw/sqrt(nuw*mag(Uw.snGrad()) + VSMALL);
    }
}

	
Foam::wallPointYPlus::yPlusCutOff = 500;
wallDistData<wallPointYPlus> yKol(mesh, ystar);

dimensionedScalar yStarLim ("ystarLim", dimLength/dimLength, 5.0);
    
volScalarField KolmogorovLimit = pos(yStarLim - yKol/ystar);
