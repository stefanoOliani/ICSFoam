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

Class
    Foam::ausmPlusUpFluxScheme

Description
    AUSM+Up flux splitting scheme as described by Liou 2006, A sequal to
    AUSM, Part II: AUSM+-up for all speeds. Included support for dynamic mesh 
    and MRF

SourceFiles
    ausmPlusUpFluxScheme.C

\*---------------------------------------------------------------------------*/

#ifndef ausmPlusUpFluxScheme_H
#define ausmPlusUpFluxScheme_H

#include "convectiveFluxScheme.H"
#include "psiThermo.H"
#include "MUSCLInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class ausmPlusUpFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class ausmPlusUpFluxScheme
:
    public convectiveFluxScheme
{
	//- Low Mach number correction
	bool lowMach_;

protected:


public:
    //- Runtime type information
    TypeName("AUSMPlusUp");


	// Constructors

        //- Construct from components
        ausmPlusUpFluxScheme
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& p
        );


    // Destructor

        ~ausmPlusUpFluxScheme();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux
		(
			surfaceScalarField& phi,
			surfaceVectorField& phiUp,
			surfaceScalarField& phiEp
		);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
