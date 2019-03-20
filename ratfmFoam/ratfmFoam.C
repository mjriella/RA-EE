/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    ratfmFoam 

Description
    A Reynolds-Averaged Two-Fluid model for turbulent fluid-particle flows.

Author
    Matthew Riella, University of Exeter.
    
References
	Riella, M., Kahraman. R., Tabor, G., (2019)
	”Inhomogeneity and anisotropy in Eulerian-Eulerian near-wall modelling“
	Int. J. Multiphase Flow, 114:9-18.

	Riella, M., Kahraman, R., and Tabor, G. R. (2018). 
	Reynolds-averaged two-fluid model prediction of moderately dilute 
	gas-solid flow over a backward-facing step. 
	International Journal of Multiphase Flow, 106:95 – 108.

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "readPPProperties.H"
    #include "initContinuityErrs.H"
	#include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "dictionary/initConvergenceCheck.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readratfmFoamControls.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaEqn.H"
            #include "dragCoeffs.H"
            #include "UEqns.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
				volScalarField nuKT = kineticTheory.nuKT();
        		#include "updateCoupling.H"
				#include "fluidTurbulence/callFluidTurbulence.H"
        		#include "updateCoupling.H"
        		#include "particleTurbulence/callParticleTurbulence.H"
            }
        }

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
