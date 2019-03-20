// convert dim scalar nu to volScalarField
// Manipulate boundaries
tmp<volScalarField> nu_field1
	(
		new volScalarField
		(
    		IOobject
    		(
        		"nu_field1",
        		runTime.timeName(),
        		mesh,
        		IOobject::NO_READ,
        		IOobject::NO_WRITE
    		),
    		mesh,
    		dimensionedScalar("nu", dimensionSet(0, 2, -1, 0, 0, 0, 0), turbulence2Properties.lookup("nu"))
		)
);

volScalarField& nu_lam1 = nu_field1();



// wall unit as defined by nu2/sqrt(taufr/rho2)
volScalarField ystar1
(
    IOobject
    (
       "ystar1",
        mesh.time().constant(),
        mesh
    ),
    mesh,
    dimensionedScalar("ysta1r", dimLength, 5.0)
);

const fvPatchList& patches1 = mesh.boundary();
forAll(patches1, patchi)
{
    if (isA<wallFvPatch>(patches1[patchi]))
    {
        const fvPatchVectorField& Uw = U2.boundaryField()[patchi];
        const scalarField& nuw = nu_lam1.boundaryField()[patchi];

        ystar1.boundaryField()[patchi] =
            nuw/sqrt(nuw*mag(Uw.snGrad()) + VSMALL);
    }
}

	
Foam::wallPointYPlus::yPlusCutOff = 500;
wallDistData<wallPointYPlus> yKol1(mesh, ystar1);

dimensionedScalar yStarLim1 ("ystarLim", dimLength/dimLength, 5.0);
    
 KolmogorovLimit = pos(yStarLim1 - yKol1/ystar1);
