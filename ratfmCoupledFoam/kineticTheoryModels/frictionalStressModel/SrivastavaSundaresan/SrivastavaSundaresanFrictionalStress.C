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

#include "SrivastavaSundaresanFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SrivastavaSundaresanFrictionalStress, 0);

    addToRunTimeSelectionTable
    (
        frictionalStressModel,
        SrivastavaSundaresanFrictionalStress,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SrivastavaSundaresanFrictionalStress::SrivastavaSundaresanFrictionalStress
(
    const dictionary& dict
)
:
    frictionalStressModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SrivastavaSundaresanFrictionalStress::~SrivastavaSundaresanFrictionalStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SrivastavaSundaresanFrictionalStress::
frictionalPressure
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{

    return
        Fr*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta)
       /pow(max(alphaMax - alpha1, scalar(5.0e-2)), p);
}


Foam::tmp<Foam::volScalarField> Foam::SrivastavaSundaresanFrictionalStress::
frictionalPressurePrime
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& alphaMax,
    const dimensionedScalar& Fr,
    const dimensionedScalar& eta,
    const dimensionedScalar& p
) const
{
    return Fr*
    (
        eta*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta - 1.0)
       *(alphaMax-alpha1)
      + p*pow(max(alpha1 - alphaMinFriction, scalar(0)), eta)
    )/pow(max(alphaMax - alpha1, scalar(5.0e-2)), p + 1.0);
}


Foam::tmp<Foam::volScalarField> Foam::SrivastavaSundaresanFrictionalStress::muf
(
    const volScalarField& alpha1,
    const dimensionedScalar& alphaMax,
    const volScalarField& pf,
    const volSymmTensorField& D,
    const dimensionedScalar& alphaMinFriction,
    const dimensionedScalar& phi
) const
{

    tmp<volScalarField> tmu
    (
        new volScalarField 
        (
            IOobject
            (
                "SrivastavaSundaresan:mu",
                alpha1.mesh().time().timeName(),
                alpha1.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            alpha1.mesh(),
            dimensionedScalar("mu", dimensionSet(1, -1, -1, 0, 0), 0.0)
        )
    );

    // Looping over all values - would be 
    // more efficient on internalField
    volScalarField& muf = tmu(); 

    forAll(D, celli)
    {
        if (alpha1[celli] > alphaMinFriction.value())
        {
            muf[celli] =
                0.5*pf[celli]*sin(phi.value())
                /(
                        sqrt((1.0/3.0)*sqr(tr(D[celli])) - invariantII(D[celli]))
                      + SMALL
                 );
        }
    }

    muf.correctBoundaryConditions();

    return tmu;
}


// ************************************************************************* //
