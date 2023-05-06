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

	defineRunTimeSelectionTable(coupledMatrix::smoother, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<coupledMatrix::smoother> coupledMatrix::smoother::New
(
        const coupledMatrix& matrix,
		const dictionary& smootherDict
)
{
    word smootherTypeName;

    smootherDict.lookup("smoother") >> smootherTypeName;

    auto* ctorPtr = dictionaryConstructorTable(smootherTypeName);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown smoother type " << smootherTypeName
            << endl << endl
            << "Valid smoother types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<smoother>(ctorPtr(matrix));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledMatrix::smoother::smoother
(
    const coupledMatrix& matrix
)
:
    matrix_(matrix)
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
