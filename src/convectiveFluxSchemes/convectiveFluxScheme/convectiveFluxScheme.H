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
    Foam::fluxScheme

Description
    Base class for inviscid flux schemes

SourceFiles
    convectiveFluxScheme.C

Authors

\*---------------------------------------------------------------------------*/

#ifndef convectiveFluxScheme_H
#define convectiveFluxScheme_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "psiThermo.H"
#include "coupledMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class convectiveFluxScheme

{
    const dictionary& dict_;
    const fvMesh& mesh_;

    const psiThermo& thermo_;

    //- References to primitive variables
    const volScalarField& rho_;
    const volVectorField& U_;
    const volScalarField& p_;

    //- Whether the class should use Lax-Friedrich approximation
    //  for the Jacobian
    bool LaxFriedrichJacob_;

    //- Face velocity in case of MRF
    surfaceScalarField* MRFFaceVelocity_;

    //- Omega field for MRF calculations
    volVectorField* MRFOmega_;


protected:
    // Protected data

    surfaceScalarField pos_;
    surfaceScalarField neg_;

private:
    // Private Member Functions

        //- Disallow copy construct
        convectiveFluxScheme(const convectiveFluxScheme&);

        //- Disallow default bitwise assignment
        void operator=(const convectiveFluxScheme&);

        void boundaryJacobian
        (
            label patchi,
            tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
            tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
            tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
        ) const;

        void addMRFSource(coupledMatrix& cMatrix) const;


public:
    //- Runtime type information
    TypeName("convectiveFluxScheme");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            convectiveFluxScheme,
            dictionary,
            (
                const dictionary& dict,
                const psiThermo& thermo,
                const volScalarField& rho,
                const volVectorField& U,
                const volScalarField& p
            ),
            (dict, thermo, rho, U, p)
        );


	// Constructors

        //- Construct from components
        convectiveFluxScheme
        (
            const word& type,
            const dictionary& dict,
			const psiThermo& thermo,
			const volScalarField& rho,
			const volVectorField& U,
            const volScalarField& p
        );


    // Selectors

        //- Return a reference to the selected gas law
        static autoPtr<convectiveFluxScheme> New
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& p
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

        surfaceScalarField& MRFFaceVelocity() const
        {
        	return *MRFFaceVelocity_;
        }

        volVectorField& MRFOmega() const
        {
        	return *MRFOmega_;
        }


    // Destructor

        virtual ~convectiveFluxScheme();


    // Member Functions

        //- Convert the boundary Jacobians w.r.t. the primitive variables and
        // add to main Jacobian
        virtual void addBoundaryTerms(coupledMatrix& cMatrix) const;

        //- Compute the temporal contribution to the approximate
        // and full Jacobian diagonal matrices
        virtual void addTemporalTerms(coupledMatrix& cMatrix, const scalarField& ddtCoeff) const;

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux
		(
			surfaceScalarField& phi,
			surfaceVectorField& phiUp,
			surfaceScalarField& phiEp
		) = 0;

        void addFluxTerms(coupledMatrix& cMatrix) const;

        virtual void addDissipationJacobian(coupledMatrix& cMatrix) const;

        void createConvectiveJacobian(coupledMatrix& cMatrix, const scalarField& ddtCoeff) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
