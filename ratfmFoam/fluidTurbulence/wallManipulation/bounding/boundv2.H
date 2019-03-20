dimensionedScalar v2Min("v2Min", sqr(dimVelocity), SMALL);

const scalar minVsfv2 = min(v2).value();

if (minVsfv2 < v2Min.value())
{
   Info<< "bounding " << v2.name()
        << ", min: " << minVsfv2
        << " max: " << max(v2).value()
        << " average: " << gAverage(v2.internalField())
        << endl;

    v2.internalField() = max
    (
        max
        (
            v2.internalField(),
            fvc::average(max(v2, v2Min))().internalField()
          * pos(-v2.internalField())
        ),
        v2Min.value()
    );

    v2.boundaryField() = max(v2.boundaryField(), v2Min.value());
}

	
