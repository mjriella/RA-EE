// correct flux
surfaceScalarField SfGradp("SfGradp", pEqn.flux()/Dp_f);

// flux correctors using rhie chow interpolation
phi1 = (fvc::interpolate(U1) & mesh.Sf())
            // pressure correction
    		+ (-SfGradp + SourceSfp)/(rho1*oneOverRatio1f)*rAU1f
            // drag correction
            + thetaD1f*(phi2.prevIter() 
                    -(fvc::interpolate(U2) & mesh.Sf()))/
                            oneOverRatio1f
            // time correction
            + thetaT1f*(phi1.oldTime() 
                    -(fvc::interpolate(U1.oldTime()) & mesh.Sf()))/
                            oneOverRatio1f; 

phi2 = (fvc::interpolate(U2) & mesh.Sf()) 
            // pressure correction
    		+ (-SfGradp + SourceSfp)/(rho2*oneOverRatio2f)*rAU2f
            // drag correction
            + thetaD2f*(phi1.prevIter()
                     -(fvc::interpolate(U1) & mesh.Sf()))/
                            oneOverRatio2f
            // time correction
            + thetaT2f*(phi2.oldTime() 
                    -(fvc::interpolate(U2.oldTime()) & mesh.Sf()))/
                            oneOverRatio2f; 

phi = alpha1f*phi1 + alpha2f*phi2;
