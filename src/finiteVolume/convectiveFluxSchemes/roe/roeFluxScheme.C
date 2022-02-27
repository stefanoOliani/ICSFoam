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

#include "addToRunTimeSelectionTable.H"
#include "roeFluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(roeFluxScheme, 0);
addToRunTimeSelectionTable(convectiveFluxScheme, roeFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::roeFluxScheme::getRoeDissipation
(
	surfaceScalarField& RoeDensity,
	surfaceVectorField& RoeVelocity,
	surfaceScalarField& RoeSoundSpeed,
	const surfaceScalarField& gammaInterp,
	surfaceScalarField& uProjRoe,
	const surfaceVectorField& normal
)
{
	RoeVelocity.dimensions().clear();
	RoeDensity.dimensions().clear();
	RoeSoundSpeed.dimensions().clear();
	uProjRoe.dimensions().clear();

	// Handy fields
	surfaceScalarField a1 = gammaInterp -1;
	surfaceScalarField theta = 0.5*a1*magSqr(RoeVelocity);
	surfaceScalarField c2 = sqr(RoeSoundSpeed);
	surfaceScalarField a2 = 1/(RoeDensity*RoeSoundSpeed*sqrt(2.0));
	surfaceScalarField a3 = RoeDensity/(RoeSoundSpeed*sqrt(2.0));
	surfaceScalarField a4 = (theta + c2)/a1;
	surfaceScalarField a5 = 1 - theta/c2;
	surfaceScalarField a6 = theta/a1;

	// Assemble inverseP
	surfaceVectorField invP11 = (normal * a5) - (RoeVelocity ^ normal)/RoeDensity;
	surfaceTensorField invP12 = a1/c2 * normal * RoeVelocity;
	invP12.component(1).ref() += normal.component(2)()/RoeDensity;
	invP12.component(2).ref() -= normal.component(1)()/RoeDensity;
	invP12.component(3).ref() -= normal.component(2)()/RoeDensity;
	invP12.component(5).ref() += normal.component(0)()/RoeDensity;
	invP12.component(6).ref() += normal.component(1)()/RoeDensity;
	invP12.component(7).ref() -= normal.component(0)()/RoeDensity;
	surfaceVectorField invP13 = -a1/c2 * normal;
	surfaceScalarField invP21 = a2 * (theta - RoeSoundSpeed*uProjRoe);
	surfaceVectorField invP22 = -a2*(a1*RoeVelocity - RoeSoundSpeed*normal);
	surfaceScalarField invP23 = a1*a2;

	surfaceScalarField invP31 = a2 * (theta + RoeSoundSpeed*uProjRoe);
	surfaceVectorField invP32 = -a2*(a1*RoeVelocity + RoeSoundSpeed*normal);
	surfaceScalarField invP33 = a1*a2;

	// Assemble |Lambda|
	surfaceScalarField Lambda1 = mag(uProjRoe);
	surfaceScalarField Lambda2 = mag(uProjRoe + RoeSoundSpeed);
	surfaceScalarField Lambda3 = mag(uProjRoe - RoeSoundSpeed);

	//Harten's entropy correction
	surfaceScalarField epsilon =  entropyFixCoeff_ * max(Lambda2,Lambda3);

    forAll(epsilon, facei)
    {
        if(Lambda1[facei] < epsilon[facei])
        {
            Lambda1[facei] = (sqr(Lambda1[facei]) + sqr(epsilon[facei]))/(2.0*epsilon[facei]);
        }

        if(Lambda2[facei] < epsilon[facei])
        {
        	Lambda2[facei] = (sqr(Lambda2[facei]) + sqr(epsilon[facei]))/(2.0*epsilon[facei]);
        }

        if(Lambda3[facei] < epsilon[facei])
        {
        	Lambda3[facei] = (sqr(Lambda3[facei]) + sqr(epsilon[facei]))/(2.0*epsilon[facei]);
        }
    }

    forAll(epsilon.boundaryField(), patchi)
    {
        const scalarField& pEpsilon = epsilon.boundaryField()[patchi];

        scalarField& pLambda1 = Lambda1.boundaryFieldRef()[patchi];
        scalarField& pLambda2 = Lambda2.boundaryFieldRef()[patchi];
        scalarField& pLambda3 = Lambda3.boundaryFieldRef()[patchi];

        forAll(pEpsilon, facei)
        {
        	 if(pLambda1[facei] < pEpsilon[facei])
        	 {
        		 pLambda1[facei] =
        				 (sqr(pLambda1[facei]) + sqr(pEpsilon[facei]))/(2.0*pEpsilon[facei]);
        	 }

        	 if(pLambda2[facei] < pEpsilon[facei])
        	 {
        		 pLambda2[facei] =
        				 (sqr(pLambda2[facei]) + sqr(pEpsilon[facei]))/(2.0*pEpsilon[facei]);
        	 }

        	 if(pLambda3[facei] < pEpsilon[facei])
        	 {
        		 pLambda3[facei] =
        				 (sqr(pLambda3[facei]) + sqr(pEpsilon[facei]))/(2.0*pEpsilon[facei]);
        	 }
        }

    }

	// Assemble P
	surfaceVectorField P11 = normal;
	surfaceScalarField P12 = a3;
	surfaceScalarField P13 = a3;
	surfaceTensorField P21 = RoeVelocity * normal;
	P21.component(1).ref() -= normal.component(2)()*RoeDensity;
	P21.component(2).ref() += normal.component(1)()*RoeDensity;
	P21.component(3).ref() += normal.component(2)()*RoeDensity;
	P21.component(5).ref() -= normal.component(0)()*RoeDensity;
	P21.component(6).ref() -= normal.component(1)()*RoeDensity;
	P21.component(7).ref() += normal.component(0)()*RoeDensity;
	surfaceVectorField P22 = a3*(RoeVelocity + RoeSoundSpeed*normal);
	surfaceVectorField P23 = a3*(RoeVelocity - RoeSoundSpeed*normal);
	surfaceVectorField P31 = normal*a6 + RoeDensity*(RoeVelocity ^ normal);
	surfaceScalarField P32 = a3*(a4 + RoeSoundSpeed*uProjRoe);
	surfaceScalarField P33 = a3*(a4 - RoeSoundSpeed*uProjRoe);

	// Compute  P x |Lambda| x inverse P and assign to corresponding fields
	dissContByRho_.set(new surfaceScalarField
			(
				((P11*Lambda1) & invP11)
				+ (Lambda2*P12*invP21)
				+ (Lambda3*P13*invP31))
			);

	dissContByRhoU_.set(new surfaceVectorField
			(
				((P11*Lambda1) & invP12)
				+ (Lambda2*P12*invP22)
				+ (Lambda3*P13*invP32))
			);
	dissContByRhoE_.set(new surfaceScalarField
				(
					((P11*Lambda1) & invP13)
					+ (Lambda2*P12*invP23)
					+ (Lambda3*P13*invP33))
				);


	dissMomByRho_.set(new surfaceVectorField
			(
				((P21*Lambda1) & invP11)
				+ (Lambda2*P22*invP21)
				+ (Lambda3*P23*invP31))
			);

	dissMomByRhoU_.set(new surfaceTensorField
			(
				((P21*Lambda1) & invP12)
				+ (Lambda2*P22*invP22)
				+ (Lambda3*P23*invP32))
			);
	dissMomByRhoE_.set(new surfaceVectorField
				(
					((P21*Lambda1) & invP13)
					+ (Lambda2*P22*invP23)
					+ (Lambda3*P23*invP33))
				);


	dissEnergyByRho_.set(new surfaceScalarField
			(
				((P31*Lambda1) & invP11)
				+ (Lambda2*P32*invP21)
				+ (Lambda3*P33*invP31))
			);

	dissEnergyByRhoU_.set(new surfaceVectorField
			(
				((P31*Lambda1) & invP12)
				+ (Lambda2*P32*invP22)
				+ (Lambda3*P33*invP32))
			);
	dissEnergyByRhoE_.set(new surfaceScalarField
				(
					((P31*Lambda1) & invP13)
					+ (Lambda2*P32*invP23)
					+ (Lambda3*P33*invP33))
				);

	dissContByRho_().dimensions().reset(dimVelocity);
	dissContByRho_().setOriented(true);
	dissContByRhoU_().dimensions().reset(dimless);
	dissContByRhoU_().setOriented(true);
	dissContByRhoE_().dimensions().reset(dimVelocity * dimMass/dimEnergy);
	dissContByRhoE_().setOriented(true);

	dissMomByRho_().dimensions().reset(dimVelocity*dimVelocity);
	dissMomByRho_().setOriented(true);
	dissMomByRhoU_().dimensions().reset(dimVelocity);
	dissMomByRhoU_().setOriented(true);
	dissMomByRhoE_().setOriented(true);

	dissEnergyByRho_().dimensions().reset(dimVelocity*dimVelocity*dimVelocity);
	dissEnergyByRho_().setOriented(true);
	dissEnergyByRhoU_().dimensions().reset(dimVelocity*dimVelocity);
	dissEnergyByRhoU_().setOriented(true);
	dissEnergyByRhoE_().dimensions().reset(dimVelocity);
	dissEnergyByRhoE_().setOriented(true);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

roeFluxScheme::roeFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p
)
:
    convectiveFluxScheme(typeName, dict, thermo, rho, U, p),
	entropyFixCoeff_(dict.subDict("convectiveFluxScheme").lookupOrDefault("entropyFixCoeff", 0.05)),
	dissContByRho_(nullptr),
	dissContByRhoU_(nullptr),
	dissContByRhoE_(nullptr),
	dissMomByRho_(nullptr),
	dissMomByRhoU_(nullptr),
	dissMomByRhoE_(nullptr),
	dissEnergyByRho_(nullptr),
	dissEnergyByRhoU_(nullptr),
	dissEnergyByRhoE_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

roeFluxScheme::~roeFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::roeFluxScheme::calcFlux
(
	surfaceScalarField& phi,
	surfaceVectorField& phiUp,
	surfaceScalarField& phiEp
)
{
    surfaceVectorField n = mesh().Sf()/mesh().magSf();

    // Prevent oriented/unoriented incompatibility below
    n.setOriented(false);

    // Left and right states
    tmp< surfaceScalarField > rho_l (fvc::interpolate(rho(), pos_, "reconstruct(rho)"));
    tmp< surfaceScalarField > rho_r (fvc::interpolate(rho(), neg_, "reconstruct(rho)"));

    tmp< surfaceScalarField > p_l = fvc::interpolate(p(), pos_, "reconstruct(rho)");
    tmp< surfaceScalarField > p_r = fvc::interpolate(p(), neg_, "reconstruct(rho)");

    surfaceVectorField U_l (IOobject("U_l", mesh().time().timeName(), mesh()), mesh(), dimensionedVector("zero", U().dimensions(), vector::zero));
    surfaceVectorField U_r (IOobject("U_r", mesh().time().timeName(), mesh()), mesh(), dimensionedVector("zero", U().dimensions(), vector::zero));
    U_l.replace(0, fvc::interpolate(U().component(0), pos_, "reconstruct(U)"));
    U_l.replace(1, fvc::interpolate(U().component(1), pos_, "reconstruct(U)"));
    U_l.replace(2, fvc::interpolate(U().component(2), pos_, "reconstruct(U)"));
    U_r.replace(0, fvc::interpolate(U().component(0), neg_, "reconstruct(U)"));
    U_r.replace(1, fvc::interpolate(U().component(1), neg_, "reconstruct(U)"));
    U_r.replace(2, fvc::interpolate(U().component(2), neg_, "reconstruct(U)"));

    // Acoustic velocity - c = sqrt(\gamma R T)
    tmp< volScalarField > gamma = thermo().gamma();

    const volScalarField& T(thermo().T());
    tmp< volScalarField > E = thermo().he(p(), T)+0.5*magSqr(U());
    tmp< surfaceScalarField > E_l (fvc::interpolate(E(), pos_, "reconstruct(T)"));
    tmp< surfaceScalarField > E_r (fvc::interpolate(E(), neg_, "reconstruct(T)"));

    // NOTE: Literature suggest enthalpy should be interpolated separately and
    // not be assembled using left and right states of energy and pressure

    tmp< volScalarField > H
    (
        (max(E(),dimensionedScalar("0", E().dimensions(), SMALL)) +
         max(p()/rho(),dimensionedScalar("0", p().dimensions()/rho().dimensions(), SMALL)))
    );
    H->rename("H");

    tmp<surfaceScalarField > H_l (fvc::interpolate(H(), pos_, "reconstruct(T)"));
    tmp<surfaceScalarField > H_r (fvc::interpolate(H(), neg_, "reconstruct(T)"));

    E.clear();
    H.clear();

    // Roe averages
    dimensionedScalar rho0("rho0",dimDensity,VSMALL);
    tmp< surfaceScalarField > coefR = sqrt(max(rho0,rho_r())/max(rho0,rho_l()));  // CHECK: Is this form of coefR correct?

    surfaceScalarField RoeDensity = coefR()*rho_l();
    surfaceVectorField RoeVelocity = (coefR()*U_r + U_l)/(coefR() + 1.0);
    tmp< surfaceScalarField > RoeEnthalpy = (coefR()*H_r() + H_l())/(coefR() + 1.0);
    tmp<surfaceScalarField> gammaInterp = fvc::interpolate(gamma());

    surfaceScalarField RoeSoundSpeed = sqrt(mag((gammaInterp()-1.0)*(RoeEnthalpy() - 0.5*magSqr(RoeVelocity))));

    gamma.clear();
    RoeEnthalpy.clear();
    coefR.clear();

    // Contravariant velocity
    tmp< surfaceScalarField > uMag_l = U_l & n;
    tmp< surfaceScalarField > uMag_r = U_r & n;
    surfaceScalarField  uProjRoe = RoeVelocity & n;

    //Conservative variables reconstruction at interface
	surfaceVectorField rhoU_l(rho_l() * U_l);
	surfaceVectorField rhoU_r(rho_r() * U_r);

	surfaceScalarField rhoE_l(rho_l() * E_l);
	surfaceScalarField rhoE_r(rho_r() * E_r);

    if (mesh().moving())
    {
    	surfaceScalarField meshVelocity = fvc::meshPhi(U())/mesh().magSf();

        meshVelocity.setOriented(false);

    	uProjRoe -= meshVelocity;
    }

    MRFFaceVelocity().setOriented(false);
    uProjRoe -= MRFFaceVelocity();

    getRoeDissipation
	(
		RoeDensity,
		RoeVelocity,
		RoeSoundSpeed,
		gammaInterp,
		uProjRoe,
		n
	);

    // Differences of conservative variables at interface
    surfaceScalarField diffRho = rho_r() - rho_l();
    surfaceVectorField diffRhoU = rhoU_r - rhoU_l;
    surfaceScalarField diffRhoE = rhoE_r - rhoE_l;

    phi = -0.5*mesh().magSf()*(dissContByRho_()*diffRho + (dissContByRhoU_() & diffRhoU) + dissContByRhoE_()*diffRhoE);
    phiUp = -0.5*mesh().magSf()*(dissMomByRho_()*diffRho + (dissMomByRhoU_() & diffRhoU) + dissMomByRhoE_()*diffRhoE);
    phiEp = -0.5*mesh().magSf()*(dissEnergyByRho_()*diffRho + (dissEnergyByRhoU_() & diffRhoU) + dissEnergyByRhoE_()*diffRhoE);

    surfaceScalarField rhoUNorm_l = rho_l()*uMag_l();
    surfaceScalarField rhoUNorm_r = rho_r()*uMag_r();

    rhoUNorm_l.setOriented(true);
    rhoUNorm_r.setOriented(true);
    n.setOriented(true);

	phi += 0.5*mesh().magSf()*(rhoUNorm_l + rhoUNorm_r);
	phiUp += 0.5*mesh().magSf()*(rhoUNorm_l*U_l + rhoUNorm_r*U_r + n*(p_l() + p_r()));
	phiEp += 0.5*mesh().magSf()*(rhoUNorm_l*H_l() +rhoUNorm_r*H_r());

    if (mesh().moving())
    {
    	surfaceScalarField meshVelocity = fvc::meshPhi(U())/mesh().magSf();

    	phi -= 0.5*mesh().magSf()*meshVelocity*(rho_l() + rho_r());
    	phiUp -= 0.5*mesh().magSf()*meshVelocity*(rhoU_l + rhoU_r);
    	phiEp -= 0.5*mesh().magSf()*meshVelocity*(rhoE_l + rhoE_r);
    }

    MRFFaceVelocity().setOriented(true);

	phi -= 0.5*mesh().magSf()*MRFFaceVelocity()*(rho_l() + rho_r());
	phiUp -= 0.5*mesh().magSf()*MRFFaceVelocity()*(rhoU_l + rhoU_r);
	phiEp -= 0.5*mesh().magSf()*MRFFaceVelocity()*(rhoE_l + rhoE_r);
}

//void Foam::roeFluxScheme::addDissipationJacobian
//(
//	coupledMatrix& cMatrix
//) const
//{
//    blockFvMatrix<scalar, scalar>& dContByRho = cMatrix.dSByS(0,0);
//    blockFvMatrix<scalar, vector>& dContByRhoU = cMatrix.dSByV(0,0);
//    blockFvMatrix<scalar, scalar>& dContByRhoE = cMatrix.dSByS(0,1);
//
//    blockFvMatrix<vector, vector>& dMomByRho = cMatrix.dVByS(0,0);
//    blockFvMatrix<vector, tensor>& dMomByRhoU = cMatrix.dVByV(0,0);
//    blockFvMatrix<vector, vector>& dMomByRhoE = cMatrix.dVByS(0,1);
//
//    blockFvMatrix<scalar, scalar>& dEnergyByRho = cMatrix.dSByS(1,0);
//    blockFvMatrix<scalar, vector>& dEnergyByRhoU = cMatrix.dSByV(1,0);
//    blockFvMatrix<scalar, scalar>& dEnergyByRhoE = cMatrix.dSByS(1,1);
//
//    dContByRho.insertDissipationBlock(dissContByRho_());
//    dContByRhoU.insertDissipationBlock(dissContByRhoU_());
//    dContByRhoE.insertDissipationBlock(dissContByRhoE_());
//
//    dMomByRho.insertDissipationBlock(dissMomByRho_());
//    dMomByRhoU.insertDissipationBlock(dissMomByRhoU_());
//    dMomByRhoE.insertDissipationBlock(dissMomByRhoE_());
//
//    dEnergyByRho.insertDissipationBlock(dissEnergyByRho_());
//    dEnergyByRhoU.insertDissipationBlock(dissEnergyByRhoU_());
//    dEnergyByRhoE.insertDissipationBlock(dissEnergyByRhoE_());
//}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
