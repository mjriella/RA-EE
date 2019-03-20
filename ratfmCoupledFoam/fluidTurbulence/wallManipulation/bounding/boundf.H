dimensionedScalar fMin("fMin", f.dimensions(), 0.0);

const scalar minVsff = min(f).value();

if (minVsff < fMin.value())
{
   Info<< "bounding " << f.name()
        << ", min: " << minVsff
        << " max: " << max(f).value()
        << " average: " << gAverage(f.internalField())
        << endl;

    f.internalField() = max
    (
        max
        (
            f.internalField(),
            fvc::average(max(f, fMin))().internalField()
          * pos(-f.internalField())
        ),
        fMin.value()
    );

    f.boundaryField() = max(f.boundaryField(), fMin.value());
}

	
