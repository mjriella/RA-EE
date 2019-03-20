/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "JohnsonJacksonParticleThetaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
//#include "twoPhaseSystem.H"
//
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        JohnsonJacksonParticleThetaFvPatchScalarField
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_(p.size()),
    specularityCoefficient_(p.size())
{}


Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const JohnsonJacksonParticleThetaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_)
{
}


Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    restitutionCoefficient_
    (
        "restitutionCoefficient",
        dimless,
        dict.lookup("restitutionCoefficient")
    ),
    specularityCoefficient_
    (
        "specularityCoefficient",
        dimless,
        dict.lookup("specularityCoefficient")
    )
{
    if
    (
        (restitutionCoefficient_.value() < 0)
     || (restitutionCoefficient_.value() > 1)
    )
    {
        FatalErrorIn
        (
            "Foam::JohnsonJacksonParticleThetaFvPatchScalarField::"
            "JohnsonJacksonParticleThetaFvPatchScalarField"
            "("
                "const fvPatch& p,"
                "const DimensionedField<scalar, volMesh>& iF,"
                "const dictionary& dict"
            ")"
        )   << "The restitution coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    if
    (
        (specularityCoefficient_.value() < 0)
     || (specularityCoefficient_.value() > 1)
    )
    {
        FatalErrorIn
        (
            "Foam::JohnsonJacksonParticleThetaFvPatchScalarField::"
            "JohnsonJacksonParticleThetaFvPatchScalarField"
            "("
                "const fvPatch& p,"
                "const DimensionedField<scalar, volMesh>& iF,"
                "const dictionary& dict"
            ")"
        )   << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const JohnsonJacksonParticleThetaFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const JohnsonJacksonParticleThetaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleThetaFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::JohnsonJacksonParticleThetaFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonParticleThetaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    
    // Amateur way - leave this redo another day 13/8/17
	// done 13/9/17

    const fvPatchScalarField& alpha =
        patch().lookupPatchField<volScalarField, scalar>("alpha1");

    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, vector>("U1");

    const fvPatchScalarField& gs0 =
    	patch().lookupPatchField<volScalarField, scalar>("gs0");

    const fvPatchScalarField& kappa =
        patch().lookupPatchField<volScalarField, scalar>("kappa");

    const fvPatchScalarField& Theta =
        patch().lookupPatchField<volScalarField, scalar>("Theta");

    const dictionary& kineticTheoryProperties = db().lookupObject<IOdictionary>
        (
            "kineticTheoryProperties"
        );
    const dimensionedScalar alphaMax(kineticTheoryProperties.lookup("alphaMax"));

    // calculate the reference value and the value fraction
    if (restitutionCoefficient_.value() != 1.0)
    {
        this->refValue() =
            (2.0/3.0)
           *specularityCoefficient_.value()
           *magSqr(U)
           /(scalar(1) - sqr(restitutionCoefficient_.value()));

        this->refGrad() = 0.0;

        scalarField c
        (
             constant::mathematical::pi
            *alpha
            *gs0
            *(scalar(1) - sqr(restitutionCoefficient_.value()))
            *sqrt(3.0*Theta)
            /max(4.0*kappa*alphaMax.value(), SMALL)
        );

        this->valueFraction() = c/(c + patch().deltaCoeffs());
    }

    // for a restitution coefficient of 1, the boundary degenerates to a fixed
    // gradient condition
    else
    {
        this->refValue() = 0.0;

        this->refGrad() =
            pos(alpha - SMALL)
           *constant::mathematical::pi
           *specularityCoefficient_.value()
           *alpha
           *gs0
           *sqrt(3.0*Theta)
           *magSqr(U)
           /max(6.0*kappa*alphaMax.value(), SMALL);


        this->valueFraction() = 0.0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleThetaFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("restitutionCoefficient")
        << restitutionCoefficient_ << token::END_STATEMENT << nl;
    os.writeKeyword("specularityCoefficient")
        << specularityCoefficient_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// ************************************************************************* //
