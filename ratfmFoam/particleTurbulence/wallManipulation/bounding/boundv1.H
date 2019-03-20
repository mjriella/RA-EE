dimensionedScalar v1Min("v1Min", sqr(dimVelocity), SMALL);

const scalar minVsfv2 = min(v1).value();

if (minVsfv2 < v1Min.value())
{
   Info<< "bounding " << v1.name()
        << ", min: " << minVsfv2
        << " max: " << max(v1).value()
        << " average: " << gAverage(v1.internalField())
        << endl;

    v1.internalField() = max
    (
        max
        (
            v1.internalField(),
            fvc::average(max(v1, v1Min))().internalField()
          * pos(-v1.internalField())
        ),
        v1Min.value()
    );

    v1.boundaryField() = max(v1.boundaryField(), v1Min.value());
}

	
