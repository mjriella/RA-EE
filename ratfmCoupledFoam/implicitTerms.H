// implicit divergence. Correction is v important. 
BlockLduSystem<vector, scalar> 
		U1corr(fvm::UDiv(correctionU1,U1,"div(U1)"));
BlockLduSystem<vector, scalar> 
		U2corr(fvm::UDiv(correctionU2,U2,"div(U2)"));

U1U2pEqn.insertBlockCoupling(6, 0, U1corr, false);
U1U2pEqn.insertBlockCoupling(6, 3, U2corr, false);

// implicit grad p. By hand as density is required.
Field<tensor7>& d = U1U2pEqn.diag().asSquare();
Field<tensor7>& u = U1U2pEqn.upper().asSquare();
Field<tensor7>& l = U1U2pEqn.lower().asSquare();
Field<vector7>& s = U1U2pEqn.source();

// Pressure weights needed for implicit treatment
tmp<surfaceInterpolationScheme<scalar> >
    tinterpScheme
    (
        surfaceInterpolationScheme<scalar>::New
        (
            p.mesh(),
            p.mesh().schemesDict().interpolationScheme("grad(p)")
        )
    );

// used for implicit linear interpolation.
const surfaceVectorField& Sf = mesh.Sf();
const unallocLabelList& owner = mesh.owner();
const unallocLabelList& neighbour = mesh.neighbour();

 // Used for sizing.
tmp<scalarField> tdiagU1 = U1Eqn.D();
scalarField& diagU1 = tdiagU1();
scalarField& upperU1 = U1Eqn.upper();
scalarField& lowerU1 = U1Eqn.lower();
vectorField& sourceU1 = U1Eqn.source();

// Pressure gradient contributions.
tmp<vectorField> tpUv = tmp<vectorField>
    (
        new vectorField(upperU1.size(), pTraits<vector>::zero)
    );
vectorField& pUv = tpUv();
tmp<vectorField> tpLv = tmp<vectorField>
    (
        new vectorField(lowerU1.size(), pTraits<vector>::zero)
    );
vectorField& pLv = tpLv();
    
tmp<vectorField> tpSv = tmp<vectorField>
    (
        new vectorField(sourceU1.size(), pTraits<vector>::zero)
    );
vectorField& pSv = tpSv();
    
tmp<vectorField> tpDv = tmp<vectorField>
    (
        new vectorField(diagU1.size(), pTraits<vector>::zero)
    );
vectorField& pDv = tpDv();
    
// Interpolation weights for use of discretisation.
tmp<surfaceScalarField> tweights = tinterpScheme().weights(p);
const surfaceScalarField& weights = tweights();

// Discretisation of grad p. 
for(int i=0; i<owner.size();i++)
{
    int o = owner[i];
    int n = neighbour[i];
    scalar w = weights.internalField()[i];
    vector s = Sf[i];
    
    pDv[o] += s*w;
    pDv[n] -= s*(1-w);
    pLv[i]  = -s*w;
    pUv[i]  = s*(1-w);
}
    	
// Get boundary condition contributions for pressure grad p. 
p.boundaryField().updateCoeffs();
forAll(p.boundaryField(),patchI)
{
        // Present fvPatchField.
        fvPatchField<scalar> & fv = p.boundaryField()[patchI];
        const fvPatch& patch = fv.patch();
        const fvsPatchScalarField& pw = weights.boundaryField()[patchI];

        // Contributions from the boundary coefficients. 
        tmp<Field<scalar> > tic = fv.valueInternalCoeffs(pw);
        Field<scalar>& ic = tic();
        tmp<Field<scalar> > tbc = fv.valueBoundaryCoeffs(pw);
        Field<scalar>& bc = tbc();
                  
        // Surface normals for this patch.
        tmp<Field<vector> > tsn = patch.Sf();
        Field<vector> sn = tsn();
    
        // Manually add the contributions from the boundary.
        forAll(fv,facei)
        {
            label c = patch.faceCells()[facei];
                       
            pDv[c] += ic[facei]*sn[facei];
            pSv[c] -= bc[facei]*sn[facei];
        }
}

// for diag
forAll(d, i)
{
    d[i](0,6) += pDv[i].x()/rho1.value();
    d[i](1,6) += pDv[i].y()/rho1.value();
    d[i](2,6) += pDv[i].z()/rho1.value();
    d[i](3,6) += pDv[i].x()/rho2.value();
    d[i](4,6) += pDv[i].y()/rho2.value();
    d[i](5,6) += pDv[i].z()/rho2.value();
}
    
// for lower
forAll(l, i)
{
    l[i](0,6) += pLv[i].x()/rho1.value();
    l[i](1,6) += pLv[i].y()/rho1.value();
    l[i](2,6) += pLv[i].z()/rho1.value();
    l[i](3,6) += pLv[i].x()/rho2.value();
    l[i](4,6) += pLv[i].y()/rho2.value();
    l[i](5,6) += pLv[i].z()/rho2.value();
}
    
// for upper
forAll(u, i)
{
    u[i](0,6) += pUv[i].x()/rho1.value();
    u[i](1,6) += pUv[i].y()/rho1.value();
    u[i](2,6) += pUv[i].z()/rho1.value();
    u[i](3,6) += pUv[i].x()/rho2.value();
    u[i](4,6) += pUv[i].y()/rho2.value();
    u[i](5,6) += pUv[i].z()/rho2.value();
}

// for source
forAll(s, i)
{
    s[i](0) += pSv[i].x()/rho1.value();
    s[i](1) += pSv[i].y()/rho1.value();
    s[i](2) += pSv[i].z()/rho1.value();
    s[i](3) += pSv[i].x()/rho2.value();
    s[i](4) += pSv[i].y()/rho2.value();
    s[i](5) += pSv[i].z()/rho2.value();
}
