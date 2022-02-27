/*---------------------------------------------------------------------------*\

    ICSFoam: a library for Implicit Coupled Simulations in OpenFOAM
  
    Copyright (C) 2022  Stefano Oliani

    https://turbofe.it

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


Author
    Stefano Oliani
    Fluid Machinery Research Group, University of Ferrara, Italy
\*---------------------------------------------------------------------------*/

#include "MRFCoupledZoneList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::MRFCoupledZoneList::zeroFilter
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphi
) const
{
    if (size())
    {
        tmp<surfaceScalarField> zphi
        (
            New
            (
                tphi,
                "zeroFilter(" + tphi().name() + ')',
                tphi().dimensions(),
                true
            )
        );

        forAll(*this, i)
        {
            operator[](i).zero(zphi.ref());
        }

        return zphi;
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


// ************************************************************************* //
