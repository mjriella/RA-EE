dimensionedScalar v2Min("v2Min", sqr(dimVelocity), SMALL);

const scalar minVsfv2 = min(v2p).value();

if (minVsfv2 < v2Min.value())
{
   Info<< "bounding " << v2p.name()
        << ", min: " << minVsfv2
        << " max: " << max(v2p).value()
        << " average: " << gAverage(v2p.internalField())
        << endl;

    v2p.internalField() = max
    (
        max
        (
            v2p.internalField(),
            fvc::average(max(v2p, v2Min))().internalField()
          * pos(-v2p.internalField())
        ),
        v2Min.value()
    );

    v2p.boundaryField() = max(v2p.boundaryField(), v2Min.value());
}

	
