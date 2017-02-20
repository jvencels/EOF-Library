/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "TableVisc.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(TableVisc, 0);
    addToRunTimeSelectionTable(viscosityModel, TableVisc, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::TableVisc::TableVisc
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(name, dimViscosity, 1)
    ),
    nuFile("constant/nuT"),
    nuGraph("nuTinterp","Temp","nuT",nuFile)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::TableVisc::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);
    correct();

    return true;
}

void Foam::viscosityModels::TableVisc::correct()
{
    //nu_ = calcNu();
    const volScalarField& T = U_.mesh().lookupObject<volScalarField>("T");

    // Iterpolate internal field
    nu_.field() = interpolateXY(T.field(), nuGraph.x(), nuGraph.y());

    /*// Iterpolate boundary field
    forAll(nu_.boundaryField(), patchi)
    {
	    nu_.boundaryField()[patchi] = interpolateXY
        (
            T.boundaryField()[patchi],
            nuGraph.x(),
            nuGraph.y()
        );
    }*/
}


// ************************************************************************* //
