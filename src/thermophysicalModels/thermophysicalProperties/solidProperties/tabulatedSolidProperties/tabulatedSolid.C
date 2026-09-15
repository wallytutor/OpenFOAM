/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "tabulatedSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tabulatedSolid, 0);
    addToRunTimeSelectionTable(solidProperties, tabulatedSolid, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Base class is initialized with meaningless coefficients, as it does not
// provide an void constructor. In the empty constructor this is the trivial
// solution, and in the dictionary constructor, this simply overrides all
// the implementation of the base class.

Foam::tabulatedSolid::tabulatedSolid()
:
    solidProperties(0, 0, 0, 0, 0),
    rhoTable_(nullptr),
    CpTable_(nullptr),
    kappaTable_(nullptr),
    hf_(0),
    emissivity_(0)
{
}

Foam::tabulatedSolid::tabulatedSolid(const dictionary& dict)
:
    solidProperties(0, 0, 0, 0, 0),
    rhoTable_
    (
        Function1<scalar>::New
        (
            "rho",
            dimTemperature,
            dimDensity,
            dict
        )
    ),
    CpTable_
    (
        Function1<scalar>::New
        (
            "Cp",
            dimTemperature,
            dimSpecificHeatCapacity,
            dict
        )
    ),
    kappaTable_
    (
        Function1<scalar>::New
        (
            dict.found("kappa") ? "kappa" : "K",
            dimTemperature,
            dimThermalConductivity,
            dict
        )
    ),
    hf_(dict.lookupBackwardsCompatible<scalar>({"Hf", "hf"})),
    emissivity_(dict.lookup<scalar>("emissivity"))
{
    solidProperties::operator=
    (
        solidProperties
        (
            rho(Tstd),
            Cp(Tstd),
            kappa(Tstd),
            hf_,
            emissivity_
        )
    );
}

Foam::tabulatedSolid::tabulatedSolid(const tabulatedSolid& ts)
:
    solidProperties(ts),
    rhoTable_(ts.rhoTable_, false),
    CpTable_(ts.CpTable_, false),
    kappaTable_(ts.kappaTable_, false),
    hf_(ts.hf_),
    emissivity_(ts.emissivity_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tabulatedSolid::write(Ostream& os) const
{
    if (rhoTable_.valid())
    {
        writeEntry(os, rhoTable_());
    }
    if (CpTable_.valid())
    {
        writeEntry(os, CpTable_());
    }
    if (kappaTable_.valid())
    {
        writeEntry(os, kappaTable_());
    }

    writeEntry(os, "Hf", hf_);
    writeEntry(os, "emissivity", emissivity_);
}

// ************************************************************************* //