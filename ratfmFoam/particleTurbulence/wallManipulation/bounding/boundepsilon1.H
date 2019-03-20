dimensionedScalar epsilon1Min("epsilon1Min", k1.dimensions()/dimTime, SMALL);    

scalar minVsf = min(epsilon1).value();

if (minVsf < epsilon1Min.value())
{
    Info<< "bounding " << epsilon1.name()
        << ", min: " << minVsf
        << " max: " << max(epsilon1).value()
        << " average: " << gAverage(epsilon1.internalField())
        << endl;

    epsilon1.internalField() = max
    (
        max
        (
            epsilon1.internalField(),
            fvc::average(max(epsilon1, epsilon1Min))().internalField()
          * pos(-epsilon1.internalField())
        ),
        epsilon1Min.value()
    );

    epsilon1.boundaryField() = max(epsilon1.boundaryField(), epsilon1Min.value());
    }

	
