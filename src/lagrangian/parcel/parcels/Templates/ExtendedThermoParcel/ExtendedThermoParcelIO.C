/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "ExtendedThermoParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::ExtendedThermoParcel<ParcelType>::propertyList_ =
    Foam::ExtendedThermoParcel<ParcelType>::propertyList();

template<class ParcelType>
const std::size_t Foam::ExtendedThermoParcel<ParcelType>::sizeofFields_
(
    sizeof(ExtendedThermoParcel<ParcelType>)
  - offsetof(ExtendedThermoParcel<ParcelType>, T_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ExtendedThermoParcel<ParcelType>::ExtendedThermoParcel(Istream& is, bool readFields)
:
    ParcelType(is, readFields),
    T_(0.0),
    Cp_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            T_ = readScalar(is);
            Cp_ = readScalar(is);
        }
        else
        {
            is.read(reinterpret_cast<char*>(&T_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check
    (
        "ExtendedThermoParcel::ExtendedThermoParcel(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::ExtendedThermoParcel<ParcelType>::readFields(CloudType& c)
{
    bool valid = c.size();

    ParcelType::readFields(c);

    IOField<scalar> T(c.fieldIOobject("T", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, T);

    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::MUST_READ), valid);
    c.checkFieldIOobject(c, Cp);


    label i = 0;
    forAllIter(typename CloudType, c, iter)
    {
        ExtendedThermoParcel<ParcelType>& p = iter();

        p.T_ = T[i];
        p.Cp_ = Cp[i];
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::ExtendedThermoParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np = c.size();

    IOField<scalar> T(c.fieldIOobject("T", IOobject::NO_READ), np);
    IOField<scalar> Cp(c.fieldIOobject("Cp", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename CloudType, c, iter)
    {
        const ExtendedThermoParcel<ParcelType>& p = iter();

        T[i] = p.T_;
        Cp[i] = p.Cp_;
        i++;
    }

    T.write(np > 0);
    Cp.write(np > 0);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::ExtendedThermoParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ExtendedThermoParcel<ParcelType>::writeFields(c);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ExtendedThermoParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.T()
            << token::SPACE << p.Cp();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.T_),
            ExtendedThermoParcel<ParcelType>::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const ExtendedThermoParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
