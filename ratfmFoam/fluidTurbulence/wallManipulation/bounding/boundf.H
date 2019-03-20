dimensionedScalar fMin("fMin", f2.dimensions(), 0.0);

const scalar minVsff = min(f2).value();

if (minVsff < fMin.value())
{
   Info<< "bounding " << f2.name()
        << ", min: " << minVsff
        << " max: " << max(f2).value()
        << " average: " << gAverage(f2.internalField())
        << endl;

    f2.internalField() = max
    (
        max
        (
            f2.internalField(),
            fvc::average(max(f2, fMin))().internalField()
          * pos(-f2.internalField())
        ),
        fMin.value()
    );

    f2.boundaryField() = max(f2.boundaryField(), fMin.value());
}

	
