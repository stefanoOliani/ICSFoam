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

#include "hllcFluxScheme.H"

#include "addToRunTimeSelectionTable.H"
#include "MUSCLInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(hllcFluxScheme, 0);
addToRunTimeSelectionTable(convectiveFluxScheme, hllcFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hllcFluxScheme::hllcFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p
)
:
    convectiveFluxScheme(typeName, dict, thermo, rho, U, p)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

hllcFluxScheme::~hllcFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::hllcFluxScheme::calcFlux
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
    tmp< volScalarField > psi = thermo().psi();
    dimensionedScalar c0("c0",dimVelocity,VSMALL);
    tmp< volScalarField > c = max(sqrt(gamma()/psi()),c0);
    tmp< surfaceScalarField > c_l = fvc::interpolate(c(), pos_, "reconstruct(T)");
    tmp< surfaceScalarField > c_r = fvc::interpolate(c(), neg_, "reconstruct(T)");

    c.clear();
    psi.clear();

    const volScalarField& T(thermo().T());
    tmp< volScalarField > E = thermo().he(p(), T)+0.5*magSqr(U());
    tmp< surfaceScalarField > E_l (fvc::interpolate(E(), pos_, "reconstruct(T)"));
    tmp< surfaceScalarField > E_r (fvc::interpolate(E(), neg_, "reconstruct(T)"));

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
    tmp< surfaceVectorField > uAvg = (coefR()*U_r + U_l)/(coefR() + 1.0);
    tmp< surfaceScalarField > HAvg = (coefR()*H_r() + H_l())/(coefR() + 1.0);
    tmp< surfaceScalarField > cAvg = sqrt(mag((fvc::interpolate(gamma())-1.0)*(HAvg() - 0.5*magSqr(uAvg()))));

    gamma.clear();
    HAvg.clear();
    coefR.clear();

    // Contravariant velocity
    tmp< surfaceScalarField > uMag_l = U_l&n;
    tmp< surfaceScalarField > uMag_r = U_r&n;
    tmp< surfaceScalarField > uMagAvg = uAvg()&n;

    if (mesh().moving())
    {
    	surfaceScalarField meshVelocity = fvc::meshPhi(U())/mesh().magSf();

        meshVelocity.setOriented(false);

        uMagAvg.ref() -= meshVelocity;

        uMag_l.ref() -= meshVelocity;
        uMag_r.ref() -= meshVelocity;

    }

    MRFFaceVelocity().setOriented(false);
    uMagAvg.ref() -= MRFFaceVelocity();

    uMag_l.ref() -= MRFFaceVelocity();
    uMag_r.ref() -= MRFFaceVelocity();

    tmp< surfaceScalarField > Sl = min(uMag_l() - c_l(), uMagAvg() - cAvg());
    tmp< surfaceScalarField > Sr = max(uMag_r() + c_r(), uMagAvg() + cAvg());

    tmp< surfaceScalarField > Sm = (rho_r()*uMag_r()*(Sr()-uMag_r()) - rho_l()*uMag_l()*(Sl()-uMag_l()) + p_l() - p_r())/
                                   (rho_r()*(Sr()-uMag_r()) - rho_l()*(Sl()-uMag_l()));
    uAvg.clear();
    cAvg.clear();
    c_l.clear();
    c_r.clear();

    // Coefficients to replace if statements
    tmp< surfaceScalarField > coefSl = pos(Sl());                           // Sl > 0 ? 1 : 0
    tmp< surfaceScalarField > coefSr = neg(Sr());                           // Sr < 0 ? 1 : 0
    tmp< surfaceScalarField > coefSm = pos(Sm());                           // Sm > 0 ? 1 : 0
    tmp< surfaceScalarField > coefSlm = (1.0 - coefSl())*coefSm();          // Sl < 0 & Sm > 0 ? 1 : 0
    tmp< surfaceScalarField > coefSmr = (1.0 - coefSm())*(1.0 - coefSr());  // Sm < 0 & Sr > 0 ? 1 : 0
    coefSm.clear();

    // Continuity
    tmp< surfaceScalarField > fluxRhoStar_l = Sm()/(Sl()-Sm())*((Sl()-uMag_l())*rho_l());
    tmp< surfaceScalarField > fluxRhoStar_r = Sm()/(Sr()-Sm())*((Sr()-uMag_r())*rho_r());

    phi =  (coefSl() *rho_l()*uMag_l() +
            coefSlm()*fluxRhoStar_l() +
            coefSmr()*fluxRhoStar_r() +
            coefSr() *rho_r()*uMag_r()
           )*mesh().magSf();
    fluxRhoStar_l.clear();
    fluxRhoStar_r.clear();

    // Momentum
    tmp< surfaceScalarField > pStar_l = rho_l()*(uMag_l()-Sl())*(uMag_l()-Sm()) + p_l();
    tmp< surfaceScalarField > pStar_r = rho_r()*(uMag_r()-Sr())*(uMag_r()-Sm()) + p_r();

    tmp< surfaceVectorField > rhoUStar_l = 1.0/(Sl()-Sm())*((Sl()-uMag_l())*rho_l()*U_l + (pStar_l() - p_l())*n);
    tmp< surfaceVectorField > rhoUStar_r = 1.0/(Sr()-Sm())*((Sr()-uMag_r())*rho_r()*U_r + (pStar_r() - p_r())*n);

    tmp< surfaceVectorField > fluxRhoUStar_l = Sm()*rhoUStar_l() + pStar_l()*n;
    tmp< surfaceVectorField > fluxRhoUStar_r = Sm()*rhoUStar_r() + pStar_r()*n;
    rhoUStar_l.clear();
    rhoUStar_r.clear();

    phiUp = (coefSl() *(rho_l()*uMag_l()*U_l + p_l()*n)+
             coefSlm()*fluxRhoUStar_l() +
             coefSmr()*fluxRhoUStar_r() +
             coefSr() *(rho_r()*uMag_r()*U_r + p_r()*n)
            )*mesh().magSf();

    fluxRhoUStar_l.clear();
    fluxRhoUStar_r.clear();

    // Energy
    tmp< surfaceScalarField > rhoEStar_l = 1.0/(Sl()-Sm())*((Sl()-uMag_l())*(rho_l()*E_l()) - p_l()*uMag_l() + pStar_l()*Sm());
    tmp< surfaceScalarField > rhoEStar_r = 1.0/(Sr()-Sm())*((Sr()-uMag_r())*(rho_r()*E_r()) - p_r()*uMag_r() + pStar_r()*Sm());
    tmp< surfaceScalarField > fluxRhoEStar_l = Sm()*(rhoEStar_l() + pStar_l()) + pStar_l() * MRFFaceVelocity();
    tmp< surfaceScalarField > fluxRhoEStar_r = Sm()*(rhoEStar_r() + pStar_r()) + pStar_r() * MRFFaceVelocity();

    if (mesh().moving())
    {
    	surfaceScalarField meshVelocity = fvc::meshPhi(U())/mesh().magSf();

        meshVelocity.setOriented(false);

    	fluxRhoEStar_l.ref() += pStar_l() * meshVelocity;
    	fluxRhoEStar_r.ref() += pStar_r() * meshVelocity;
    }

    MRFFaceVelocity().setOriented(true);

    phi.setOriented(true);

    phiEp = (coefSl() *rho_l()*H_l()*uMag_l() +
             coefSlm()*fluxRhoEStar_l() +
             coefSmr()*fluxRhoEStar_r() +
             coefSr() *rho_r()*H_r()*uMag_r()
            )*mesh().magSf();

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
