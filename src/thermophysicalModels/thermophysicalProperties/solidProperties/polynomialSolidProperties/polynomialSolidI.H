/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::scalar Foam::polynomialSolid::rho() const
{
    return rho(Tstd);
}


inline Foam::scalar Foam::polynomialSolid::rho(const scalar T) const
{
    return rhoCoeffs_.value(T);
}


inline Foam::scalar Foam::polynomialSolid::Cp() const
{
    return Cp(Tstd);
}


inline Foam::scalar Foam::polynomialSolid::Cp(const scalar T) const
{
    return CpCoeffs_.value(T);
}


inline Foam::scalar Foam::polynomialSolid::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::polynomialSolid::hs(const scalar T) const
{
    return CpCoeffs_.integral(Tstd, T);
}


inline Foam::scalar Foam::polynomialSolid::ha(const scalar T) const
{
    return hs(T) + hf();
}


inline Foam::scalar Foam::polynomialSolid::es(const scalar T) const
{
    return hs(T);
}


inline Foam::scalar Foam::polynomialSolid::ea(const scalar T) const
{
    return es(T) + hf();
}


inline Foam::scalar Foam::polynomialSolid::kappa() const
{
    return kappa(Tstd);
}


inline Foam::scalar Foam::polynomialSolid::kappa(const scalar T) const
{
    return kappaCoeffs_.value(T);
}


inline Foam::scalar Foam::polynomialSolid::emissivity() const
{
    return emissivity_;
}

// ************************************************************************* //