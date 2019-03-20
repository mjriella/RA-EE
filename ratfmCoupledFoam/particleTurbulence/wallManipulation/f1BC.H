
{
    labelList cellBoundaryFaceCount(f1.size(), 0);

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall f fields to zero
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                f1[faceCelli] = symmTensor::zero;
            }
        }
    }

    //- Accumulate the wall face contributions to f
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

         if (isType<wallFvPatch>(curPatch))
        {

            if ((f1.boundaryField()[patchi]).type() != "fixedInternalValue")
            {
                FatalErrorIn("boundary-condition evaluation")
                    << f1.boundaryField()[patchi].type()
                    << " is the wrong f1 patchField type for"
                    << " boundary-condition on patch "
                    << curPatch.name() << nl
                    << "    should be fixedInternalValue"
                    << exit(FatalError);
            }

            const scalarField& nuw = nu2.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                // For corner cells (with two boundary or more faces),
                // f in the near-wall cell is calculated
                // as an average

                cellBoundaryFaceCount[faceCelli]++;

                symmTensor fbc =
                  - 20.0*sqr(nuw[facei])*R1[faceCelli]/epsilon1[faceCelli]
                  / sqr(sqr(y[patchi][facei]));

                f1[faceCelli].yy() += fbc.yy();
                f1[faceCelli].xy() += fbc.xy();
                f1[faceCelli].yz() += fbc.yz();

                f1[faceCelli].xx() += -0.5*fbc.yy();
                f1[faceCelli].zz() += -0.5*fbc.yy();

                f1[faceCelli].xz() *= 0.0;
            }
        }
    }


    // Perform the averaging
    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                f1[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}


// ************************************************************************* //
