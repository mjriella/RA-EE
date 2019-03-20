    if (runTime.outputTime())
    {
        volVectorField Ur
        (
            IOobject
            (
                "Ur",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U1 - U2
        );

        volScalarField Ka = K;

        runTime.write();
    }

	// Useful for debug.        
	/*Info << "Max velocity, U1 " << max(U1).value() << endl;
    Info << "Mean velocity, U1 " << average(U1).value() << endl;
    Info << "Max velocity, U2 " << max(U2).value() << endl;
    Info << "Mean velocity, U2 " << average(U2).value() << endl;
    Info << "Pressure " << max(p).value() << endl;*/
