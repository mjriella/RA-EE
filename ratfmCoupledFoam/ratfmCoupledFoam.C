/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    ratfmCoupledFoam 

Description
    Fully-coupled Reynolds-Averaged Two-Fluid model for the solution of 
    fluid-particle systems.

Author
    Matthew Riella, University of Exeter.
    
References
 	Riella, M., Kahraman, R., Tabor, G. (2019).
 	Fully-coupled pressure-based two-fluid solver for the solution of 
	turbulent fluid-particle systems.
	Computer and Fluids.

 	Riella, M. (2019).
	Turbulence modelling of fluid-particle interaction. 
	University of Exeter. PhD Thesis.

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "fixedValueFvsPatchFields.H"
#include "nearWallDist.H"
#include "wallDistData.H"
#include "wallPointYPlus.H"
#include "wallFvPatch.H"
#include "dragModel.H"
#include "phaseModel.H"
#include "kineticTheoryModel.H"
#include "pimpleControl.H"
#include "fvBlockMatrix.H"
#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "readPPProperties.H"
#   include "createTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"
#   include "initContinuityErrs.H"
#   include "dictionary/initConvergenceCheck.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "dictionary/readratfmControls.H"
#       include "dictionary/readBlockSolverControls.H"
#       include "dictionary/readFieldBounds.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        U1.storePrevIter();
        U2.storePrevIter();
        phi1.storePrevIter();
        phi2.storePrevIter();
        p.storePrevIter();

        // Initialize the U1U2p block.
        fvBlockMatrix<vector7> U1U2pEqn(U1U2p);

        // solve volume fraction eq.
#       include "alphaEqn.H"
        // construct drag.
#		include "dragCalc.H"
        // discretise and insert mom eqs.
#       include "U1U2.H"
        // construct flux corrs and pressure.
# 		include "cFluxPressure.H"
        // remove drag in source and make implicit.
#		include "couplingCoeffs.H"
		// implicit div (with drag corr) and grad p.
#		include "implicitTerms.H"

        // Solve the block matrix.
        maxResidual = cmptMax(U1U2pEqn.solve().initialResidual());
    
        // Retrieve solution.
    	U1U2pEqn.retrieveSolution(0, U1.internalField());
    	U1U2pEqn.retrieveSolution(3, U2.internalField());
    	U1U2pEqn.retrieveSolution(6, p.internalField());

        // Correct updated boundaries here.
        U1.correctBoundaryConditions();
        U2.correctBoundaryConditions();
        p.correctBoundaryConditions();

        // Recalculate mixture velocity.
        U = alpha1*U1 + alpha2*U2; 
        U.correctBoundaryConditions();

		// Update fluxes after solution.
#		include "updateFluxes.H"

		// Evaluate cont errors and bound if necs.
#		include "continuityErrs.H"
#       include "dictionary/boundPU.H"

        p.relax();

        // Solve phase energy equations.
		volScalarField nuKT = kineticTheory.nuKT();
        #include "updateCoupling.H"
        #include "particleTurbulence/callParticleTurbulence.H"
		#include "fluidTurbulence/callFluidTurbulence.H"

        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "dictionary/convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
