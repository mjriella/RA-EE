    IOdictionary ppProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    scalar preAlphaExp(readScalar(ppProperties.lookup("preAlphaExp")));

    scalar alphaMax(readScalar(ppProperties.lookup("alphaMax")));

    scalar expMax(readScalar(ppProperties.lookup("expMax")));

    dimensionedScalar g0(ppProperties.lookup("g0"));
