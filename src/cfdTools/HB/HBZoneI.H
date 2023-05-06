/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2022 Stefano Oliani
-------------------------------------------------------------------------------
License
    This file is part of ICSFOAM.

    ICSFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

inline const Foam::word& Foam::HBZone::name() const
{
    return name_;
}


inline bool Foam::HBZone::active() const
{
    return active_;
}


inline const Foam::scalarList& Foam::HBZone::omegaList() const
{
    return omegaList_;
}

inline const Foam::vector& Foam::HBZone::rotationAxis() const
{
    return rotationAxis_;
}

inline const Foam::point& Foam::HBZone::rotationCentre() const
{
    return rotationCentre_;
}

inline const Foam::RectangularMatrix<complex>& Foam::HBZone::E() const
{
    return E_;
}

inline const Foam::RectangularMatrix<complex>& Foam::HBZone::EInv() const
{
    return EInv_;
}

inline const Foam::scalarSquareMatrix& Foam::HBZone::D() const
{
    return D_;
}


// ************************************************************************* //
