dimensionedScalar f1Min("f1Min", f1.dimensions(), 0.0);

const scalar minVsff = min(f1).value();

if (minVsff < f1Min.value())
{
   Info<< "bounding " << f1.name()
        << ", min: " << minVsff
        << " max: " << max(f1).value()
        << " average: " << gAverage(f1.internalField())
        << endl;

    f1.internalField() = max
    (
        max
        (
            f1.internalField(),
            fvc::average(max(f1, f1Min))().internalField()
          * pos(-f1.internalField())
        ),
        f1Min.value()
    );

    f1.boundaryField() = max(f1.boundaryField(), f1Min.value());
}

	
