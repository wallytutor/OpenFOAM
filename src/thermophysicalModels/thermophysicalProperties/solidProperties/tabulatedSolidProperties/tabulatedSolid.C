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
    kappaTable_(nullptr)
{
}

Foam::tabulatedSolid::tabulatedSolid(const dictionary& dict)
:
    solidProperties(0, 0, 0, 0, 0),
    rhoTable_("rho", dimless, dict.subDict("rho")),
    CpTable_("Cp", dimless, dict.subDict("Cp")),
    kappaTable_("kappa", dimless, dict.subDict("kappa")),
    hf_(dict.lookup<scalar>("Hf")),
    emissivity_(dict.lookup<scalar>("emissivity"))
{}

// tabulatedSolid::tabulatedSolid(const tabulatedSolid& ms)
// :
//     solidProperties(0, 0, 0, 0, 0),
//     rhoTable_(ms.rhoTable_.valid() ? ms.rhoTable_->clone() : nullptr),
//     CpTable_(ms.CpTable_.valid() ? ms.CpTable_->clone() : nullptr),
//     kappaTable_(ms.kappaTable_.valid() ? ms.kappaTable_->clone() : nullptr),
//     hf_(ms.hf_),
//     emissivity_(ms.emissivity_),
// {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tabulatedSolid::write(Ostream& os) const
{
    // if (rhoTable_.valid())
    // {
    //     rhoTable_->write(os);
    // }
    // if (CpTable_.valid())
    // {
    //     CpTable_->write(os);
    // }
    // if (kappaTable_.valid())
    // {
    //     kappaTable_->write(os);
    // }

    os << hf_ << token::SPACE << emissivity_;
}

// ************************************************************************* //