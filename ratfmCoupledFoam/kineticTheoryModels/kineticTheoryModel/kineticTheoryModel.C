/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kineticTheoryModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const Foam::phaseModel& phase1,
    const Foam::volVectorField& U2,
    const Foam::volScalarField& alpha1,
    const Foam::dragModel& drag1
)
:
    phase1_(phase1),
    U1_(phase1.U()),
    U2_(U2),
    alpha1_(alpha1),
    phi1_(phase1.phi()),
    drag1_(drag1),
    rho1_(phase1.rho()),
    da_(phase1.d()),
    nu1_(phase1.nu()),
    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            U1_.time().constant(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    kineticTheory_(kineticTheoryProperties_.lookup("kineticTheory")),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    conductivityModel_
    (
        conductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            kineticTheoryProperties_
        )
    ),
    granularPressureModel_
    (
        granularPressureModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_(kineticTheoryProperties_.lookup("e")),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    Theta_
    (
        IOobject
        (
            "Theta",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh()
    ),
    mu1_
    (
        IOobject
        (
            "mu1",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    nuKT_
    (
        IOobject
        (
            "nuKT",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),
    pa_
    (
        IOobject
        (
            "pa",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    gs0_
    (
        IOobject
        (
            "gs0",
            U1_.time().timeName(),
            U1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U1_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModel::solve()
{
    if (!kineticTheory_)
    {
        return;
    }

    //- coupling lookup
    volScalarField epsilon1 = U1_.mesh().lookupObject<volScalarField> ("epsilon1");
    dimensionedScalar alphaSmall("alphaSmall", dimless, 1.0e-6);

    // - Declare local refs for KT equation
    surfaceScalarField phi(rho1_*phi1_*fvc::interpolate(alpha1_));

    volTensorField dU = fvc::grad(U1_);
    volTensorField dUT = dU.T();  
    volTensorField D = 0.5*(dU + dUT);
   
    // - particle turb visocisty needed in conductivity 
    volScalarField nut1 = U1_.mesh().lookupObject<volScalarField> ("nut1");
    volScalarField mut = nut1*rho1_;

    // - Turbulent prandtl number - 0.85
    volScalarField mut1 = 1.79*mut*alpha1_;

    const scalar sqrtPi = sqrt(mathematicalConstant::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));
 
    // NB, drag has alpha in equation here
    volScalarField Ur(mag(U1_ - U2_));
    volScalarField alpha2Prim(alpha1_*(1.0 - alpha1_)*drag1_.K(Ur));
	
	// radial dist
    gs0_ = radialModel_->g0
    (
        min(max(alpha1_, scalar(1e-6)), alphaMax_ - 0.01),
        alphaMax_
    );

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha1_,
            gs0_,
            rho1_,
            e_
        )
    );

    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha1_, Theta_, gs0_, rho1_, da_, e_);

    // particle viscosity (Table 3.2, p.47)
    mu1_ = viscosityModel_->mu1(alpha1_, Theta_, gs0_, rho1_, da_, e_);

    dimensionedScalar Tsmall
    (
        "small",
        dimensionSet(0 , 2 ,-2 ,0 , 0, 0, 0),
        1.0e-6
    );

    volScalarField ThetaSqrt(sqrt(Theta_));

    // dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
    (
        12.0*(1.0 - sqr(e_))*sqr(alpha1_)*rho1_*gs0_*(1.0/da_)*ThetaSqrt/sqrtPi
    );

    // Coupling from drag
    volScalarField J1 = 3.0*alpha2Prim;

    // Build stress tensor like fluid
    // compP is the compression term which acts similiar to bulk
    tmp<volTensorField> tgradU = fvc::grad(U1_);
    volScalarField tau = 2.0*mu1_*(tgradU() && dev(twoSymm(tgradU())));
    volScalarField compP = PsCoeff*fvc::div(U1_);
    tgradU.clear();

    // Granular temperature transport equation
    fvScalarMatrix ThetaEqn
    (
        1.5*
        (
        	fvm::ddt(alpha1_*rho1_, Theta_)
      	  + fvm::div(phi, Theta_, "div(phi1,Theta)")
        )
      - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
      - fvm::laplacian(mut1, Theta_, "laplacian(mut1,Theta)")
      ==
        tau
      - fvm::SuSp(compP, Theta_)
      + fvm::Sp(-gammaCoeff, Theta_)
      + fvm::Sp(-J1, Theta_)  
      + alpha1_*rho1_*epsilon1
    );
    ThetaEqn.relax();
    ThetaEqn.solve();

    // Limit for stability 
    Theta_.max(1.0e-15);
    Theta_.min(1e6); 

    // Update after solution of theta
    PsCoeff = granularPressureModel_->granularPressureCoeff
        (
            alpha1_,
            gs0_,
            rho1_,
            e_
        );

    mu1_ = viscosityModel_->mu1(alpha1_, Theta_, gs0_, rho1_, da_, e_);

    // Frictional stress has been removed as it shouldnt be here
    
    // update particle pressure
    pa_ = PsCoeff*Theta_;
    mu1_.min(1.0e+2);
    mu1_.max(0.0);

    // Update viscosity
    nuKT_ = mu1_/rho1_;
}



// ************************************************************************* /}/
