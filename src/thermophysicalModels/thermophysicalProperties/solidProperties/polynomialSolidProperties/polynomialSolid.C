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

#include "polynomialSolid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polynomialSolid, 0);
    addToRunTimeSelectionTable(solidProperties, polynomialSolid, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Base class is initialized with meaningless coefficients, as it does not
// provide an void constructor. In the empty constructor this is the trivial
// solution, and in the dictionary constructor, this simply overrides all
// the implementation of the base class.

Foam::polynomialSolid::polynomialSolid()
:
    solidProperties(0, 0, 0, 0, 0),
    rhoCoeffs_(),
    CpCoeffs_(),
    kappaCoeffs_(),
    hf_(0),
    emissivity_(0)
{
}

Foam::polynomialSolid::polynomialSolid(const dictionary& dict)
:
    solidProperties(0, 0, 0, 0, 0),
    rhoCoeffs_(dict.lookup("rhoCoeffs")),
    CpCoeffs_(dict.lookup("CpCoeffs")),
    kappaCoeffs_(dict.lookup(dict.found("kappaCoeffs") ? "kappaCoeffs" : "KCoeffs")),
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

Foam::polynomialSolid::polynomialSolid(const polynomialSolid& ps)
:
    solidProperties(ps),
    rhoCoeffs_(ps.rhoCoeffs_),
    CpCoeffs_(ps.CpCoeffs_),
    kappaCoeffs_(ps.kappaCoeffs_),
    hf_(ps.hf_),
    emissivity_(ps.emissivity_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polynomialSolid::write(Ostream& os) const
{
    writeEntry(os, "rhoCoeffs", rhoCoeffs_);
    writeEntry(os, "CpCoeffs", CpCoeffs_);
    writeEntry(os, "kappaCoeffs", kappaCoeffs_);
    writeEntry(os, "Hf", hf_);
    writeEntry(os, "emissivity", emissivity_);
}

// ************************************************************************* //