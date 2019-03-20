    
dimensionedScalar k1Min("k1Min", sqr(dimVelocity), SMALL);

const scalar minVsfk1 = min(k1).value();

if (minVsfk1 < k1Min.value())
{
   Info<< "bounding " << k1.name()
        << ", min: " << minVsfk1
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

	
