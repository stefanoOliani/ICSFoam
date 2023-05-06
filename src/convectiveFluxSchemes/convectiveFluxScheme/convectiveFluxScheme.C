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
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "blockFvOperators.H"
#include "addToRunTimeSelectionTable.H"

#include "MUSCLInterpolation.H"
#include "MUSCLInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(convectiveFluxScheme, 0);
defineRunTimeSelectionTable(convectiveFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void convectiveFluxScheme::boundaryJacobian
(
    label patchi,
    tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
    tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
    tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
) const
{
    const volScalarField& T(thermo_.T());
    tmp< volScalarField > tcv = thermo_.Cv();
    const volScalarField& cv = tcv();

    volScalarField rhoEn = rho_*(thermo_.he(p_, T)+0.5*magSqr(U_));

    const volScalarField::Boundary& pbf = p_.boundaryField();
    const volVectorField::Boundary& ubf = U_.boundaryField();
    const volScalarField::Boundary& tbf = T.boundaryField();

    const scalarField& wP = mesh_.weights().boundaryField()[patchi];

    // ValueInternalCoeffs
    // fixedValue   -> 0
    // zeroGradient -> 1

    tmp < vectorField > tuVIC = ubf[patchi].valueInternalCoeffs(wP);
    const vectorField& uVIC = tuVIC();

    tmp < scalarField > tpVIC = pbf[patchi].valueInternalCoeffs(wP);
    const scalarField& pVIC = tpVIC();

    tmp < scalarField > ttVIC = tbf[patchi].valueInternalCoeffs(wP);
    const scalarField& tVIC = ttVIC();

    const vectorField& SfB = mesh_.Sf().boundaryField()[patchi];

    const scalarField& rhoB = rho_.boundaryField()[patchi];

    tmp<vectorField> UB = ubf[patchi];

    tmp<scalarField> UrelBdotSf = UB() & SfB;

    if (mesh_.moving())
    {
        const surfaceScalarField& meshPhi(mesh_.phi());
        UrelBdotSf.ref() -= meshPhi.boundaryField()[patchi];
    }

    UrelBdotSf.ref() -= MRFFaceVelocity_->boundaryField()[patchi]*mesh_.magSf().boundaryField()[patchi];

    const scalarField& pB = pbf[patchi];
    const scalarField& TB = tbf[patchi];

    const scalarField& rhoEB = rhoEn.boundaryField()[patchi];

    const scalarField& cvB = cv.boundaryField()[patchi];

    dContFluxdp = rhoB/pB * UrelBdotSf() * pVIC;
    dContFluxdU = rhoB*cmptMultiply(SfB,uVIC);
    dContFluxdT = -rhoB/TB * UrelBdotSf() * tVIC;

    dMomFluxdp = (rhoB/pB*UB() * UrelBdotSf() + SfB) * pVIC;
    dMomFluxdU = rhoB * UB() * cmptMultiply(SfB,uVIC);
    tmp<vectorField> dMomFluxdUDiag = rhoB * UrelBdotSf() * uVIC;
    dMomFluxdT = -rhoB/TB*UB() * UrelBdotSf() * tVIC;

    dEnergyFluxdp = (rhoEB/pB*UrelBdotSf() + (UB() & SfB)) * pVIC;
    dEnergyFluxdU = cmptMultiply(SfB, uVIC)*(rhoEB+pB) + rhoB*UrelBdotSf()*cmptMultiply(UB(), uVIC);
    dEnergyFluxdT = UrelBdotSf*(rhoB*cvB-rhoEB/TB)*tVIC;

    dMomFluxdU.ref().replace(tensor::XX, dMomFluxdU().component(tensor::XX)+dMomFluxdUDiag().component(vector::X));
    dMomFluxdU.ref().replace(tensor::YY, dMomFluxdU().component(tensor::YY)+dMomFluxdUDiag().component(vector::Y));
    dMomFluxdU.ref().replace(tensor::ZZ, dMomFluxdU().component(tensor::ZZ)+dMomFluxdUDiag().component(vector::Z));

}


void convectiveFluxScheme::addMRFSource(coupledMatrix& cMatrix) const
{
	volVectorField rhoU = rho_ * U_;

	cMatrix.dVByV(0,0).source() -=
			(MRFOmega().primitiveField() ^ rhoU.primitiveField()) * mesh_.V();

	tensorField& MomDiag = cMatrix.dVByV(0,0).diag();

	MomDiag.component(1).ref() -= MRFOmega().primitiveField().component(2)() * mesh_.V();
	MomDiag.component(2).ref() += MRFOmega().primitiveField().component(1)() * mesh_.V();
	MomDiag.component(3).ref() += MRFOmega().primitiveField().component(2)() * mesh_.V();
	MomDiag.component(5).ref() -= MRFOmega().primitiveField().component(0)() * mesh_.V();
	MomDiag.component(6).ref() -= MRFOmega().primitiveField().component(1)() * mesh_.V();
	MomDiag.component(7).ref() += MRFOmega().primitiveField().component(0)() * mesh_.V();

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

convectiveFluxScheme::convectiveFluxScheme
(
    const word& type,
    const dictionary& dict,
	const psiThermo& thermo,
	const volScalarField& rho,
	const volVectorField& U,
    const volScalarField& p
)
:
	dict_(dict),
	mesh_(U.mesh()),
	thermo_(thermo),
	rho_(rho),
	U_(U),
	p_(p),
	LaxFriedrichJacob_(true),
	MRFFaceVelocity_(nullptr),
	MRFOmega_(nullptr),
    pos_(surfaceScalarField
    (
        IOobject
        (
            "pos",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("pos", dimless, 1.0)
    )),
    neg_(surfaceScalarField
    (
        IOobject
        (
            "neg",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("neg", dimless, -1.0)
    ))
{
	MRFFaceVelocity_ = new surfaceScalarField
						(
							IOobject
							(
								"MRFFaceVelocity",
								mesh_.time().timeName(),
								mesh_
							),
							mesh_,
							dimensionedScalar("zero", dimVelocity, Zero)
						);

	MRFOmega_ = new volVectorField
		(
			IOobject
			(
				"MRFOmega",
				mesh_.time().timeName(),
				mesh_
			),
			mesh_,
			dimensionedVector("zero", dimless/dimTime, Zero)
		);
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

convectiveFluxScheme::~convectiveFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void convectiveFluxScheme::addBoundaryTerms(coupledMatrix& cMatrix) const
{
    tmp< volScalarField > tcv = thermo_.Cv();
    const volScalarField& cv = tcv();
    tmp< volScalarField > tgamma(thermo_.gamma());
    const volScalarField& gamma = tgamma();

    const volScalarField& T(thermo_.T());

    blockFvMatrix<scalar, scalar>& dContByRho = cMatrix.dSByS(0,0);
    blockFvMatrix<scalar, vector>& dContByRhoU = cMatrix.dSByV(0,0);
    blockFvMatrix<scalar, scalar>& dContByRhoE = cMatrix.dSByS(0,1);
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

            volScalarField rhoEn = rho_*(thermo_.he(p_, T)+0.5*magSqr(U_));

            tmp<scalarField> rhoEI = rhoEn.boundaryField()[patchi].patchInternalField();
            tmp<scalarField> gammaI = gamma.boundaryField()[patchi].patchInternalField();
            tmp<scalarField> cvI = cv.boundaryField()[patchi].patchInternalField();

            tmp<scalarField> dPdRho = 0.5*(gammaI()-1)*magSqr(UI());
            tmp<vectorField> dUdRho = -UI()/rhoI();
            tmp<scalarField> dTdRho = -1.0/(cvI()*rhoI())*(rhoEI()/rhoI()-magSqr(UI()));

            tmp<vectorField> dPdRhoU = -(gammaI()-1)*UI();
            tmp<sphericalTensorField> dUdRhoU = 1.0/rhoI()*sphericalTensor::I;
            tmp<vectorField> dTdRhoU = -UI()/(cvI()*rhoI());

            tmp<scalarField> dPdRhoE = gammaI()-1;
            tmp<vectorField> dUdRhoE(new vectorField(mesh_.boundary()[patchi].size(), vector::zero));
            tmp<scalarField> dTdRhoE = 1.0/(cvI()*rhoI());

            scalarField& diagContByRho = dContByRho.diag();
            vectorField& diagContByRhoU = dContByRhoU.diag();
            scalarField& diagContByRhoE = dContByRhoE.diag();

            vectorField& diagMomByRho = dMomByRho.diag();
            tensorField& diagMomByRhoU = dMomByRhoU.diag();
            vectorField& diagMomByRhoE = dMomByRhoE.diag();

            scalarField& diagEnergyByRho = dEnergyByRho.diag();
            vectorField& diagEnergyByRhoU = dEnergyByRhoU.diag();
            scalarField& diagEnergyByRhoE = dEnergyByRhoE.diag();

            forAll(patch, bfacei)
            {

                label iIntCell = patch.faceCells()[bfacei];

                // Continuity coupled to rho
                diagContByRho[iIntCell] += dContFluxdp()[bfacei] * dPdRho()[bfacei];
                diagContByRho[iIntCell] += dContFluxdU()[bfacei] & dUdRho()[bfacei];
                diagContByRho[iIntCell] += dContFluxdT()[bfacei] * dTdRho()[bfacei];

                // Continuity coupled to rhoU
                diagContByRhoU[iIntCell] += dContFluxdp()[bfacei] * dPdRhoU()[bfacei];
                diagContByRhoU[iIntCell] += dContFluxdU()[bfacei] & dUdRhoU()[bfacei];
                diagContByRhoU[iIntCell] += dContFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Continuity coupled to rhoE
                diagContByRhoE[iIntCell] += dContFluxdp()[bfacei] * dPdRhoE()[bfacei];
                diagContByRhoE[iIntCell] += dContFluxdU()[bfacei] & dUdRhoE()[bfacei];
                diagContByRhoE[iIntCell] += dContFluxdT()[bfacei] * dTdRhoE()[bfacei];


                // Momentum coupled to rho
                diagMomByRho[iIntCell] += dMomFluxdp()[bfacei] * dPdRho()[bfacei];
                diagMomByRho[iIntCell] += dMomFluxdU()[bfacei] & dUdRho()[bfacei];
                diagMomByRho[iIntCell] += dMomFluxdT()[bfacei] * dTdRho()[bfacei];

                // Momentum coupled to rhoU
                diagMomByRhoU[iIntCell] += dMomFluxdp()[bfacei] * dPdRhoU()[bfacei];
                diagMomByRhoU[iIntCell] += dMomFluxdU()[bfacei] & dUdRhoU()[bfacei];
                diagMomByRhoU[iIntCell] += dMomFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Momentum coupled to rhoE
                diagMomByRhoE[iIntCell] += dMomFluxdp()[bfacei] * dPdRhoE()[bfacei];
                diagMomByRhoE[iIntCell] += dMomFluxdU()[bfacei] & dUdRhoE()[bfacei];
                diagMomByRhoE[iIntCell] += dMomFluxdT()[bfacei] * dTdRhoE()[bfacei];


                // Energy coupled to rho
                diagEnergyByRho[iIntCell] += dEnergyFluxdp()[bfacei] * dPdRho()[bfacei];
                diagEnergyByRho[iIntCell] += dEnergyFluxdU()[bfacei] & dUdRho()[bfacei];
                diagEnergyByRho[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRho()[bfacei];

                // Energy coupled to rhoU
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdp()[bfacei] * dPdRhoU()[bfacei];
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdU()[bfacei] & dUdRhoU()[bfacei];
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Energy coupled to rhoE
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdp()[bfacei] * dPdRhoE()[bfacei];
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdU()[bfacei] & dUdRhoE()[bfacei];
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRhoE()[bfacei];

            }
        }
    }
}


void convectiveFluxScheme::addTemporalTerms(coupledMatrix& cMatrix, const scalarField& ddtCoeff) const
{
    blockFvMatrix<scalar, scalar>& dContByRho = cMatrix.dSByS(0,0);
    blockFvMatrix<vector, tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
    blockFvMatrix<scalar, scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);

    // Temporal contribution (Weak form)
    scalarField diagCoeff = ddtCoeff*mesh_.V();
    dContByRho.diag() += diagCoeff;
    dMomByRhoU.diag() += diagCoeff*tensor::I;
    dEnergyByRhoE.diag() += diagCoeff;
}


// Assemble full matrix with diagonal and off-diagonal contribution from the
// flux terms.
void convectiveFluxScheme::addFluxTerms(coupledMatrix& cMatrix) const
{
    surfaceScalarField w(mesh_.weights());

    w.setOriented(false);  // Workaround; it should not be oriented

    tmp< volScalarField > tgamma(thermo_.gamma());
    const volScalarField& gamma = tgamma();

    // Left and right states
    surfaceVectorField U_l (IOobject("U_l", mesh().time().timeName(), mesh()), mesh(), dimensionedVector("zero", U().dimensions(), vector::zero));
    surfaceVectorField U_r (IOobject("U_r", mesh().time().timeName(), mesh()), mesh(), dimensionedVector("zero", U().dimensions(), vector::zero));

    U_l.replace(0, fvc::interpolate(U().component(0), pos_, "reconstruct(U)"));
    U_l.replace(1, fvc::interpolate(U().component(1), pos_, "reconstruct(U)"));
    U_l.replace(2, fvc::interpolate(U().component(2), pos_, "reconstruct(U)"));
    U_r.replace(0, fvc::interpolate(U().component(0), neg_, "reconstruct(U)"));
    U_r.replace(1, fvc::interpolate(U().component(1), neg_, "reconstruct(U)"));
    U_r.replace(2, fvc::interpolate(U().component(2), neg_, "reconstruct(U)"));

    surfaceScalarField gamma_interp (fvc::interpolate(gamma));

    const volScalarField& T = thermo_.T();
    volScalarField  E = thermo_.he(p_, T)+0.5*magSqr(U_);

    surfaceScalarField E_l (fvc::interpolate(E, pos_, "reconstruct(T)"));
    surfaceScalarField E_r (fvc::interpolate(E, neg_, "reconstruct(T)"));

    surfaceScalarField theta_l = 0.5*(gamma_interp - 1) * (magSqr(U_l));
    surfaceScalarField theta_r = 0.5*(gamma_interp - 1) * (magSqr(U_r));

    surfaceScalarField a1_l = gamma_interp * E_l - theta_l;
    surfaceScalarField a1_r = gamma_interp * E_r - theta_r;

    surfaceScalarField a2_l = gamma_interp - 1;
    surfaceScalarField a2_r = gamma_interp - 1;

    surfaceVectorField n = mesh().Sf()/mesh().magSf();

    // Prevent oriented/unoriented incompatibility below
    n.setOriented(false);

    surfaceScalarField projU_l = U_l & n;
    surfaceScalarField projU_r = U_r & n;

    blockFvMatrix<scalar, scalar>& dContByRho = cMatrix.dSByS(0,0);
    blockFvMatrix<scalar, vector>& dContByRhoU = cMatrix.dSByV(0,0);
    cMatrix.dSByS(0,1).diag();

    blockFvMatrix<vector, vector>& dMomByRho = cMatrix.dVByS(0,0);
    blockFvMatrix<vector, tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
    blockFvMatrix<vector, vector>& dMomByRhoE = cMatrix.dVByS(0,1);

    blockFvMatrix<scalar, scalar>& dEnergyByRho = cMatrix.dSByS(1,0);
    blockFvMatrix<scalar, vector>& dEnergyByRhoU = cMatrix.dSByV(1,0);
    blockFvMatrix<scalar, scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);

//	// Diagonal and off-diagonal contribution of convective part
	dContByRhoU.insertBlock(n,n);

	dMomByRho.insertBlock
	(
		(n*theta_l - U_l*projU_l)(),
		(n*theta_r - U_r*projU_r)()
	);

	dMomByRhoU.insertBlock
	(
		(U_l*n - a2_l*n*U_l + projU_l*tensor::I)(),
		(U_r*n - a2_r*n*U_r + projU_r*tensor::I)()
	);

	dMomByRhoE.insertBlock(n*a2_l, n*a2_r);

	dEnergyByRho.insertBlock
	(
		(projU_l*(theta_l-a1_l))(),
		(projU_r*(theta_r-a1_r))()
	);

	dEnergyByRhoU.insertBlock
	(
		(n*a1_l - a2_l*U_l*projU_l)(),
		(n*a1_r - a2_r*U_r*projU_r)()
	);

	dEnergyByRhoE.insertBlock
	(
		(gamma_interp*projU_l)(),
		(gamma_interp*projU_r)()
	);

    // Moving mesh part
    if (mesh_.moving())
    {
        const surfaceScalarField& meshPhi(mesh_.phi());
        tmp<blockFvMatrix<scalar,scalar> > divMeshPhi = fvj::div(w, meshPhi);

        dContByRho -= divMeshPhi();
        dMomByRhoU -= divMeshPhi()*tensor::I;
        dEnergyByRhoE -= divMeshPhi;
    }

    tmp<blockFvMatrix<scalar,scalar> > MRFdivMeshPhi = fvj::div(w, MRFFaceVelocity()*mesh_.magSf());

    dContByRho -= MRFdivMeshPhi();
    dMomByRhoU -= MRFdivMeshPhi()*tensor::I;
    dEnergyByRhoE -= MRFdivMeshPhi;

    addDissipationJacobian(cMatrix);
}


void convectiveFluxScheme::addDissipationJacobian(coupledMatrix& cMatrix) const
{
	// This function is called by derived flux classes if LaxFriedrichJacob_ = false

    blockFvMatrix<scalar, scalar>& dContByRho = cMatrix.dSByS(0,0);
    blockFvMatrix<vector, tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
    blockFvMatrix<scalar, scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);

    tmp<volScalarField> tGamma = thermo_.gamma();
    const volScalarField& gamma = tGamma();

    // Wave speed: Lax-Friedrich flux approximation of left-hand side Jacobian
    tmp< volScalarField > c = sqrt(gamma/thermo_.psi());

	surfaceScalarField lambdaConv
	(
		IOobject
		(
			"lambdaConv",
			mesh_.time().timeName(),
			mesh_
		),
		mesh_,
		dimensionedScalar("lambdaConv", dimVelocity, 0.0)
	);

	if (mesh_.moving())
	{
		lambdaConv = (fvc::interpolate(c)
						+ mag((fvc::interpolate(U_) & mesh_.Sf()/mesh_.magSf())
						- MRFFaceVelocity()
						- fvc::meshPhi(U_)/mesh_.magSf()));
	}
	else
	{
		lambdaConv = (fvc::interpolate(c)
						+ mag((fvc::interpolate(U_) & mesh_.Sf()/mesh_.magSf())
						- MRFFaceVelocity()
					  ));
	}

	dContByRho.insertDissipationBlock(lambdaConv);

	dMomByRhoU.insertDissipationBlock(lambdaConv *tensor::I);

	dEnergyByRhoE.insertDissipationBlock(lambdaConv);

}


void convectiveFluxScheme::createConvectiveJacobian(coupledMatrix& cMatrix, const scalarField& ddtCoeff) const
{
	addFluxTerms(cMatrix);

	addBoundaryTerms(cMatrix);

	addTemporalTerms(cMatrix, ddtCoeff);

	addMRFSource(cMatrix);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
