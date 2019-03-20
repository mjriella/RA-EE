// Need coupling coefficients here for in parallel 
volScalarField couplingInU1(drag12);
scalarField inU1 = couplingInU1.internalField();
volScalarField couplingInU2(drag21);
scalarField inU2 = couplingInU2.internalField();

// drag off diagonal coeffs
U1U2pEqn.insertEquationCoupling(0,3, -inU1);
U1U2pEqn.insertEquationCoupling(1,4, -inU1);
U1U2pEqn.insertEquationCoupling(2,5, -inU1);
U1U2pEqn.insertEquationCoupling(3,0, -inU2);
U1U2pEqn.insertEquationCoupling(4,1, -inU2);
U1U2pEqn.insertEquationCoupling(5,2, -inU2);

// remove drag in mom
U1U2pEqn.updateSourceCoupling();
