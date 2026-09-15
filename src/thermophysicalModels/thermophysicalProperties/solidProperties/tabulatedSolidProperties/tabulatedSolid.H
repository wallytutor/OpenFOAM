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
    Foam::tabulatedSolid

Description
    Temperature dependent solid properties

SourceFiles
    tabulatedSolid.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedSolid_H
#define tabulatedSolid_H

#include "solidProperties.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class tabulatedSolid Declaration
\*---------------------------------------------------------------------------*/

class tabulatedSolid
:
    public solidProperties
{
    // Private Data
        //- Tabulated/Function1 density rho(T)
        autoPtr<Function1<scalar>> rhoTable_;

        //- Tabulated/Function1 specific heat capacity Cp(T)
        autoPtr<Function1<scalar>> CpTable_;

        //- Tabulated/Function1 thermal conductivity kappa(T)
        autoPtr<Function1<scalar>> kappaTable_;

        //- Heat of formation [J/kg]
        scalar hf_;

        //- Emissivity
        scalar emissivity_;
public:

    //- Runtime type information
    TypeName("tabulatedSolid");

    // Constructors

        //- Construct null
        tabulatedSolid();

        //- Construct from dictionary
        tabulatedSolid(const dictionary& dict);

        //- Copy constructor
        tabulatedSolid(const tabulatedSolid& ts);

        //- Construct and return a clone
        virtual autoPtr<solidProperties> clone() const
        {
            return autoPtr<solidProperties>(new tabulatedSolid(*this));
        }

    // Member Functions

        // Physical constants which define the solidProperties

            //- Density [kg/m^3]
            inline scalar rho() const;
            inline scalar rho(const scalar T) const;

            //- Specific heat capacity [J/kg/K]
            inline scalar Cp() const;
            inline scalar Cp(const scalar T) const;

            //- Heat of formation [J/kg]
            inline scalar hf() const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar T) const;

            //- Absolute enthalpy [J/kg]
            inline scalar ha(const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar es(const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar ea(const scalar T) const;

            //- Thermal conductivity [W/m/K]
            inline scalar kappa() const;
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

#include "tabulatedSolidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //