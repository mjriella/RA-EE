    
dimensionedScalar k1Min("k1Min", sqr(dimVelocity), SMALL);

const scalar minVsfk = min(k1).value();

if (minVsfk < k1Min.value())
{
   Info<< "bounding " << k1.name()
        << ", min: " << minVsfk
        << " max: " << max(k1).value()
        << " average: " << gAverage(k1.internalField())
        << endl;

    k1.internalField() = max
    (
        max
        (
            k1.internalField(),
            fvc::average(max(k1, k1Min))().internalField()
          * pos(-k1.internalField())
        ),
        k1Min.value()
    );

    k1.boundaryField() = max(k1.boundaryField(), k1Min.value());
}

	
