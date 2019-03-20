
{
    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& p = patches[patchi];

        if (isA<wallFvPatch>(p))
        {
            epsEqn.setValues
            (
                p.faceCells(),
                epsilon2.boundaryField()[patchi].patchInternalField()
            );
        }
    }
}
