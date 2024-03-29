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

#include "../phaseLagAMILduInterfaceField/phaseLagAMILduInterfaceField.H"

#include "diagTensorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(phaseLagAMILduInterfaceField, 0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseLagAMILduInterfaceField::~phaseLagAMILduInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::phaseLagAMILduInterfaceField::transformCoupleField
(
    solveScalarField& f,
    const direction cmpt
) const
{
    if (doTransform())
    {
        if (forwardT().size() == 1)
        {
            f *= pow(diag(forwardT()[0]).component(cmpt), rank());
        }
        else
        {
            f *= pow(diag(forwardT())().component(cmpt), rank());
        }
    }
}


// ************************************************************************* //
