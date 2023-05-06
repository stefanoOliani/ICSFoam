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

#include "convectiveFluxScheme.H"
#include "error.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<convectiveFluxScheme> convectiveFluxScheme::New
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p
)
{
    word fluxSchemeTypeName;

    dict.subDict("convectiveFluxScheme").lookup("fluxScheme") >> fluxSchemeTypeName;

    Info << "Selecting flux scheme " << fluxSchemeTypeName << endl;

    auto* ctorPtr = dictionaryConstructorTable(fluxSchemeTypeName);

    if (!ctorPtr)
    {
        FatalErrorIn
        (
            "convectiveFluxScheme::New()"
        )   << "Unknown fluxScheme type " << fluxSchemeTypeName
            << endl << endl
            << "Valid convectiveFluxScheme types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<convectiveFluxScheme>(ctorPtr(dict, thermo, rho, U, p));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
