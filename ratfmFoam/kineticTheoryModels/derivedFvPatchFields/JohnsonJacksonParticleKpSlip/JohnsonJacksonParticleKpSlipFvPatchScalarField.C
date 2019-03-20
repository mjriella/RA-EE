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

#include "JohnsonJacksonParticleKpSlipFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        JohnsonJacksonParticleKpSlipFvPatchScalarField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::
JohnsonJacksonParticleKpSlipFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    partialSlipFvPatchScalarField(p, iF),
    specularityCoefficient_(p.size())
{}


Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::
JohnsonJacksonParticleKpSlipFvPatchScalarField
(
    const JohnsonJacksonParticleKpSlipFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchScalarField(ptf, p, iF, mapper),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::
JohnsonJacksonParticleKpSlipFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchScalarField(p, iF),
    specularityCoefficient_
    (
        "specularityCoefficient",
        dimless,
        dict.lookup("specularityCoefficient")
    )
{
    if
    (
        (specularityCoefficient_.value() < 0)
     || (specularityCoefficient_.value() > 1)
    )
    {
        FatalErrorIn
        (
            "("
                "Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::"
                "JohnsonJacksonParticleKpSlipFvPatchScalarField"
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


Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::
JohnsonJacksonParticleKpSlipFvPatchScalarField
(
    const JohnsonJacksonParticleKpSlipFvPatchScalarField& ptf
)
:
    partialSlipFvPatchScalarField(ptf),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::
JohnsonJacksonParticleKpSlipFvPatchScalarField
(
    const JohnsonJacksonParticleKpSlipFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    partialSlipFvPatchScalarField(ptf, iF),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    partialSlipFvPatchScalarField::autoMap(m);
}


void Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    partialSlipFvPatchScalarField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    
    // Amateur way - leave this redo another day 13/8/17
    //
    // done 13/9/17

   const fvPatchScalarField& alpha =
        patch().lookupPatchField<volScalarField, scalar>("alpha1");

    const fvPatchScalarField& nu =
        patch().lookupPatchField<volScalarField, scalar>("nuKT");

    const fvPatchScalarField& gs0 =
    	patch().lookupPatchField<volScalarField, scalar>("gs0");


    const fvPatchScalarField& Theta =
        patch().lookupPatchField<volScalarField, scalar>("Theta");

    const dictionary& kineticTheoryProperties = db().lookupObject<IOdictionary>
        (
            "kineticTheoryProperties"
        );
    const dimensionedScalar alphaMax(kineticTheoryProperties.lookup("alphaMax"));

    //Info << "particle visco = kp BC " << max(nu) << endl;


    // calculate the slip value fraction
    scalarField c
    (
        constant::mathematical::pi
       *alpha
       *gs0
       *specularityCoefficient_.value()
       *sqrt(6.0*Theta)
       /max(3.0*nu*alphaMax.value(), SMALL)
    );


    this->valueFraction() = c/(c + patch().deltaCoeffs());

    partialSlipFvPatchScalarField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleKpSlipFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("specularityCoefficient")
        << specularityCoefficient_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// ************************************************************************* //
