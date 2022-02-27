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

#include "viscousFluxScheme.H"

#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "blockFvOperators.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void viscousFluxScheme::boundaryJacobian
(
    label patchi,
    tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
    tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
    tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
) const
{
    const volScalarField& T(thermo_.T());

    surfaceScalarField muEff = fvc::interpolate(turbulence_.muEff());
    surfaceScalarField alphaEff = fvc::interpolate(turbulence_.alphaEff());

    const volVectorField::Boundary& ubf = U_.boundaryField();
    const volScalarField::Boundary& tbf = T.boundaryField();

    const scalarField gammaB = thermo_.gamma()->boundaryField()[patchi];
	const scalarField cvB = thermo_.Cv()->boundaryField()[patchi];

	const vectorField uGIC = ubf[patchi].gradientInternalCoeffs();
	const scalarField TGIC = tbf[patchi].gradientInternalCoeffs();

	const scalarField& magSfB = mesh_.magSf().boundaryField()[patchi];

	dContFluxdp = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
	dContFluxdU = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
	dContFluxdT = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));

	dMomFluxdp = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));

	dMomFluxdU = tmp<tensorField>(new tensorField(mesh_.boundary()[patchi].patch().size(), Zero));

	tmp<vectorField> dMomFluxdUDiag(-muEff.boundaryField()[patchi]*uGIC*magSfB);

	dMomFluxdU.ref().replace(tensor::XX, dMomFluxdUDiag().component(vector::X));
	dMomFluxdU.ref().replace(tensor::YY, dMomFluxdUDiag().component(vector::Y));
	dMomFluxdU.ref().replace(tensor::ZZ, dMomFluxdUDiag().component(vector::Z));

	dMomFluxdT = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));

	dEnergyFluxdp = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
	dEnergyFluxdU = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));

	dEnergyFluxdT = -alphaEff.boundaryField()[patchi]*cvB*TGIC*magSfB;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

viscousFluxScheme::viscousFluxScheme
(
    const dictionary& dict,
	const psiThermo& thermo,
	const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
	const compressible::turbulenceModel& turbulence
)
:
	dict_(dict),
	mesh_(rho.mesh()),
	thermo_(thermo),
	rho_(rho),
	U_(U),
	p_(p),
	LaxFriedrichJacob_(dict.subDict("viscousFluxScheme").lookupOrDefault("LaxFriedrichJacobian", true)),
	turbulence_(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

viscousFluxScheme::~viscousFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void viscousFluxScheme::addBoundaryTerms(coupledMatrix& cMatrix) const
{
    tmp< volScalarField > tcv = thermo_.Cv();
    const volScalarField& cv = tcv();

    blockFvMatrix<vector, vector>& dMomByRho = cMatrix.dVByS(0,0);
    blockFvMatrix<vector, tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
    blockFvMatrix<vector, vector>& dMomByRhoE = cMatrix.dVByS(0,1);
    blockFvMatrix<scalar, scalar>& dEnergyByRho = cMatrix.dSByS(1,0);
    blockFvMatrix<scalar, vector>& dEnergyByRhoU = cMatrix.dSByV(1,0);
    blockFvMatrix<scalar, scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);

    // Calculate diagonal contribution of the geometric boundaries

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];

        if(!patch.coupled() && !isA<emptyFvPatch>(patch))
        {

            tmp<scalarField> dContFluxdp;
            tmp<vectorField> dContFluxdU;
            tmp<scalarField> dContFluxdT;

            tmp<vectorField> dMomFluxdp;
            tmp<tensorField> dMomFluxdU;
            tmp<vectorField> dMomFluxdT;

            tmp<scalarField> dEnergyFluxdp;
            tmp<vectorField> dEnergyFluxdU;
            tmp<scalarField> dEnergyFluxdT;

            boundaryJacobian
            (
                patchi,
                dContFluxdp, dContFluxdU, dContFluxdT,
                dMomFluxdp, dMomFluxdU, dMomFluxdT,
                dEnergyFluxdp, dEnergyFluxdU, dEnergyFluxdT
            );


            tmp<scalarField> rhoI = rho_.boundaryField()[patchi].patchInternalField();
            tmp<vectorField> UI = U_.boundaryField()[patchi].patchInternalField();

            const volScalarField& T(thermo().T());
            tmp< volScalarField > E = thermo().he(p_, T)+0.5*magSqr(U_);

            tmp<scalarField> EI = E().boundaryField()[patchi].patchInternalField();
            tmp<scalarField> cvI = cv.boundaryField()[patchi].patchInternalField();

            tmp<vectorField> dUdRho = -UI()/rhoI();
            tmp<scalarField> dTdRho = -1.0/(cvI()*rhoI())*(EI()-magSqr(UI()));

            tmp<sphericalTensorField> dUdRhoU = 1.0/rhoI()*sphericalTensor::I;
            tmp<vectorField> dTdRhoU = -UI()/(cvI()*rhoI());

            tmp<vectorField> dUdRhoE(new vectorField(mesh_.boundary()[patchi].size(), vector::zero));
            tmp<scalarField> dTdRhoE = 1.0/(cvI()*rhoI());

            vectorField& diagMomByRho = dMomByRho.diag();
            tensorField& diagMomByRhoU = dMomByRhoU.diag();
            vectorField& diagMomByRhoE = dMomByRhoE.diag();

            scalarField& diagEnergyByRho = dEnergyByRho.diag();
            vectorField& diagEnergyByRhoU = dEnergyByRhoU.diag();
            scalarField& diagEnergyByRhoE = dEnergyByRhoE.diag();

            forAll(patch, bfacei)
            {

                label iIntCell = patch.faceCells()[bfacei];

                // Momentum coupled to rho
                diagMomByRho[iIntCell] += dMomFluxdU()[bfacei] & dUdRho()[bfacei];

                // Momentum coupled to rhoU
                diagMomByRhoU[iIntCell] += dMomFluxdU()[bfacei] & dUdRhoU()[bfacei];

                // Momentum coupled to rhoE
                diagMomByRhoE[iIntCell] += dMomFluxdU()[bfacei] & dUdRhoE()[bfacei];


                // Energy coupled to rho
                diagEnergyByRho[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRho()[bfacei];

                // Energy coupled to rhoU
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Energy coupled to rhoE
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRhoE()[bfacei];

            }
        }
    }
}


// Assemble full matrix with diagonal and off-diagonal contribution from the
// flux terms.
void viscousFluxScheme::addFluxTerms(coupledMatrix& cMatrix) const
{
    surfaceScalarField muEff = fvc::interpolate(turbulence_.muEff());
    surfaceScalarField alphaEff = fvc::interpolate(turbulence_.alphaEff());

    blockFvMatrix<scalar,scalar>& dContByRho = cMatrix.dSByS(0,0);
    blockFvMatrix<vector,vector>& dMomByRho = cMatrix.dVByS(0,0);
    blockFvMatrix<vector,tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
    blockFvMatrix<scalar,scalar>& dEnergyByRho = cMatrix.dSByS(1,0);
    blockFvMatrix<scalar,vector>& dEnergyByRhoU = cMatrix.dSByV(1,0);
    blockFvMatrix<scalar,scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);

    // Add Jacobian for the Laplacian part of the viscous terms

    if (LaxFriedrichJacob_)
    {
        // Viscous wave speed (Luo et al, 2001)
        // CHECK: Blazek proposes an alternative approximation of viscous spectral
        // radius

        surfaceScalarField lambdaVisc = (muEff+alphaEff)/fvc::interpolate(rho_);

        tmp<blockFvMatrix<scalar,scalar> > stab = fvj::laplacian(0.5*lambdaVisc, geometricOneField());

        dContByRho -= stab();
        dMomByRhoU -= stab()*tensor::I;
        dEnergyByRhoE -= stab();
    }
    else
    {
    	scalar dummyVal = 0.0;

        const volScalarField& T(thermo().T());
        tmp< volScalarField > E = thermo().he(p_, T)+0.5*magSqr(U_);

        dMomByRho -= fvj::laplacian(muEff, -U_/rho_);
        dMomByRhoU -= fvj::laplacian(muEff, 1.0/rho_)*tensor::I;

        dEnergyByRho -= fvj::laplacian(alphaEff, -E()/rho_+magSqr(U_)/rho_);
        dEnergyByRhoU -= fvj::laplacian(alphaEff, -U_/rho_, dummyVal);
        dEnergyByRhoE -= fvj::laplacian(alphaEff, 1.0/rho_);
    }
}


void viscousFluxScheme::createViscousJacobian(coupledMatrix& cMatrix) const
{
	addFluxTerms(cMatrix);

	if (!LaxFriedrichJacob_)
	{
		addBoundaryTerms(cMatrix);
	}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
