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

Class
    Foam::roeFluxScheme

Description
    Roe flux splitting scheme with suport for dynamic mesh and MRF

SourceFiles
    roeFluxScheme.C

Authors


\*---------------------------------------------------------------------------*/

#ifndef roeFluxScheme_H
#define roeFluxScheme_H

#include "convectiveFluxScheme.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class roeFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class roeFluxScheme
:
    public convectiveFluxScheme
{

	//- Coefficient for Harten's entropy fix
	const scalar entropyFixCoeff_;

	//- Components of Roe Dissipation Matrix |A|
    autoPtr<surfaceScalarField> dissContByRho_;
    autoPtr<surfaceVectorField> dissContByRhoU_;
    autoPtr<surfaceScalarField>  dissContByRhoE_;

    autoPtr<surfaceVectorField> dissMomByRho_;
    autoPtr<surfaceTensorField> dissMomByRhoU_;
    autoPtr<surfaceVectorField> dissMomByRhoE_;

    autoPtr<surfaceScalarField> dissEnergyByRho_;
    autoPtr<surfaceVectorField> dissEnergyByRhoU_;
    autoPtr<surfaceScalarField> dissEnergyByRhoE_;

    void getRoeDissipation
	(
		surfaceScalarField& RoeDensity,
		surfaceVectorField& RoeVelocity,
		surfaceScalarField& RoeSoundSpeed,
		const surfaceScalarField& gammaInterp,
		surfaceScalarField& uProjRoe,
		const surfaceVectorField& normal
	);



protected:

public:
    //- Runtime type information
    TypeName("ROE");


	// Constructors

        //- Construct from components
        roeFluxScheme
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& p
        );


    // Destructor

        ~roeFluxScheme();


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
