dimensionedScalar k2Min("k2Min", sqr(dimVelocity), SMALL);

const scalar minVsfk = min(k2).value();

if (minVsfk < k2Min.value())
{
   Info<< "bounding " << k2.name()
        << ", min: " << minVsfk
        << " max: " << max(k2).value()
        << " average: " << gAverage(k2.internalField())
        << endl;

    k2.internalField() = max
    (
        max
        (
            k2.internalField(),
            fvc::average(max(k2, k2Min))().internalField()
          * pos(-k2.internalField())
        ),
        k2Min.value()
    );

    k2.boundaryField() = max(k2.boundaryField(), k2Min.value());
}

	
