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

inline Foam::scalar Foam::tabulatedSolid::rho() const
{
    return rho(Tstd);
}


inline Foam::scalar Foam::tabulatedSolid::rho(const scalar T) const
{
    return rhoTable_.valid() ? rhoTable_->value(T) : 0.0;
}


inline Foam::scalar Foam::tabulatedSolid::Cp() const
{
    return Cp(Tstd);
}


inline Foam::scalar Foam::tabulatedSolid::Cp(const scalar T) const
{
    return CpTable_.valid() ? CpTable_->value(T) : 0.0;
}


inline Foam::scalar Foam::tabulatedSolid::hf() const
{
    return hf_;
}


inline Foam::scalar Foam::tabulatedSolid::hs(const scalar T) const
{
    return CpTable_.valid() ? CpTable_->integral(Tstd, T) : 0.0;
}


inline Foam::scalar Foam::tabulatedSolid::ha(const scalar T) const
{
    return hs(T) + hf();
}


inline Foam::scalar Foam::tabulatedSolid::es(const scalar T) const
{
    return hs(T);
}


inline Foam::scalar Foam::tabulatedSolid::ea(const scalar T) const
{
    return es(T) + hf();
}


inline Foam::scalar Foam::tabulatedSolid::kappa() const
{
    return kappa(Tstd);
}


inline Foam::scalar Foam::tabulatedSolid::kappa(const scalar T) const
{
    return kappaTable_.valid() ? kappaTable_->value(T) : 0.0;
}


inline Foam::scalar Foam::tabulatedSolid::emissivity() const
{
    return emissivity_;
}

// ************************************************************************* //