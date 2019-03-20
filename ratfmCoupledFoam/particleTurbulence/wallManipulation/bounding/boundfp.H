dimensionedScalar fMin("fMin", fp.dimensions(), 0.0);

const scalar minVsff = min(f).value();

if (minVsff < fMin.value())
{
   Info<< "bounding " << fp.name()
        << ", min: " << minVsff
        << " max: " << max(fp).value()
        << " average: " << gAverage(fp.internalField())
        << endl;

    fp.internalField() = max
    (
        max
        (
            fp.internalField(),
            fvc::average(max(fp, fMin))().internalField()
          * pos(-fp.internalField())
        ),
        fMin.value()
    );

    fp.boundaryField() = max(fp.boundaryField(), fMin.value());
}

	
