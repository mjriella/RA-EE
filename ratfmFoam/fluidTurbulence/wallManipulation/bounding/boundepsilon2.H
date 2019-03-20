dimensionedScalar epsilon2Min("epsilon2Min", k2.dimensions()/dimTime, SMALL);    

scalar minVsf = min(epsilon2).value();

if (minVsf < epsilon2Min.value())
{
    Info<< "bounding " << epsilon2.name()
        << ", min: " << minVsf
        << " max: " << max(epsilon2).value()
        << " average: " << gAverage(epsilon2.internalField())
        << endl;

    epsilon2.internalField() = max
    (
        max
        (
            epsilon2.internalField(),
            fvc::average(max(epsilon2, epsilon2Min))().internalField()
          * pos(-epsilon2.internalField())
        ),
        epsilon2Min.value()
    );

    epsilon2.boundaryField() = max(epsilon2.boundaryField(), epsilon2Min.value());
    }

	
