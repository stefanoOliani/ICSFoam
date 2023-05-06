/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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

#include "error.H"

#include "coupledMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineRunTimeSelectionTable(coupledMatrix::preconditioner, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<coupledMatrix::preconditioner> coupledMatrix::preconditioner::New
(
        const solver& sol,
        const dictionary& dict
)
{
    word preconditionerTypeName;

    dict.lookup("preconditioner") >> preconditionerTypeName;

    auto* ctorPtr = dictionaryConstructorTable(preconditionerTypeName);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown preconditioner type " << preconditionerTypeName
            << endl << endl
            << "Valid preconditioner types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<preconditioner>(ctorPtr(sol, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
