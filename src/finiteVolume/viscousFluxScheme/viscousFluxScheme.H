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

#ifndef viscousFluxScheme_H
#define viscousFluxScheme_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "psiThermo.H"
#include "coupledMatrix.H"
#include "turbulentFluidThermoModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class viscousFluxScheme

{
    const dictionary& dict_;
    const fvMesh& mesh_;

    const psiThermo& thermo_;

    //- References to conservative variables
    const volScalarField& rho_;
    const volVectorField& U_;
    const volScalarField& p_;

    //- Whether the class should use Lax-Friedrich approximation
    //  for the Jacobian
    bool LaxFriedrichJacob_;

    const compressible::turbulenceModel& turbulence_;


protected:


private:
    // Private Member Functions

        //- Disallow copy construct
        viscousFluxScheme(const viscousFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const viscousFluxScheme&);

        void boundaryJacobian
        (
            label patchi,
            tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
            tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
            tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
        ) const;


public:

	// Constructors

        //- Construct from components
        viscousFluxScheme
        (
            const dictionary& dict,
			const psiThermo& thermo,
			const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& p,
			const compressible::turbulenceModel& turbulence
        );


     // Access

        const dictionary& dict() const
        {
        	return dict_;
        }

        const fvMesh& mesh() const
        {
        	return mesh_;
        }

        const psiThermo& thermo() const
        {
        	return thermo_;
        }

        const volScalarField& rho() const
        {
        	return rho_;
        }

        const volVectorField& U() const
        {
        	return U_;
        }

        const volScalarField& p() const
        {
        	return p_;
        }


    // Destructor

        virtual ~viscousFluxScheme();


    // Member Functions

        //- Convert the boundary Jacobians w.r.t. the primitive variables and
        // add to main Jacobian
        virtual void addBoundaryTerms(coupledMatrix& cMatrix) const;

        void addFluxTerms(coupledMatrix& cMatrix) const;

        void createViscousJacobian(coupledMatrix& cMatrix) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
