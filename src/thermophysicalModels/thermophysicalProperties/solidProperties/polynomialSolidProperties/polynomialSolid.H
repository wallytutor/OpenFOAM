/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::polynomialSolid

Description
    Temperature dependent solid properties

SourceFiles
    polynomialSolid.C

\*---------------------------------------------------------------------------*/

#ifndef polynomialSolid_H
#define polynomialSolid_H

#include "solidProperties.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polynomialSolid Declaration
\*---------------------------------------------------------------------------*/

class polynomialSolid
:
    public solidProperties
{
    // Private Data
        // Polynomial coefficients for rho(T) = r0 + r1*T + r2*T^2 + ...
        Polynomial<8> rhoCoeffs_;

        // Polynomial coefficients for Cp(T) = c0 + c1*T + c2*T^2 + ...
        Polynomial<8> CpCoeffs_;

        // Polynomial coefficients for kappa(T) = c0 + c1*T + c2*T^2 + ...
        Polynomial<8> kappaCoeffs_;

        //- Heat of formation [J/kg]
        scalar hf_;

        //- Emissivity
        scalar emissivity_;
public:

    //- Runtime type information
    TypeName("polynomialSolid");

    // Constructors

        //- Construct null
        polynomialSolid();

        //- Construct from dictionary
        polynomialSolid(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<solidProperties> clone() const
        {
            return autoPtr<solidProperties>(new polynomialSolid(*this));
        }

    // Member Functions

        // Physical constants which define the solidProperties

            //- Density [kg/m^3]
            inline scalar rho(const scalar T) const;

            //- Specific heat capacity [J/kg/K]
            inline scalar Cp(const scalar T) const;

            //- Heat of formation [J/kg]
            inline scalar hf() const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar T) const;

            //- Absolute enthalpy [J/kg]
            inline scalar ha(scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar es(const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar ea(scalar T) const;

            //- Thermal conductivity [W/m/K]
            inline scalar kappa(const scalar T) const;

            //- Emissivity []
            inline scalar emissivity() const;

    // I-O

        //- Read and set the properties present it the given dictionary
        // void readIfPresent(const dictionary& dict);

        //- Write the custolSolid properties
        virtual void write(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const solidProperties& s);
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polynomialSolidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
