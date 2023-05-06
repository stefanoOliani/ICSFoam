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

#include "ausmPlusUpFluxScheme.H"

#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "cellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ausmPlusUpFluxScheme, 0);
addToRunTimeSelectionTable(convectiveFluxScheme, ausmPlusUpFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ausmPlusUpFluxScheme::ausmPlusUpFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p
)
:
    convectiveFluxScheme(typeName, dict, thermo, rho, U, p),
	lowMach_(dict.subDict("convectiveFluxScheme").lookupOrDefault<Switch>("lowMachAusm", true))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

ausmPlusUpFluxScheme::~ausmPlusUpFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ausmPlusUpFluxScheme::calcFlux
(
	surfaceScalarField& phi,
	surfaceVectorField& phiUp,
	surfaceScalarField& phiEp
)
{
    surfaceVectorField U_L (IOobject("U_l", mesh().time().timeName(), mesh()), mesh(), dimensionedVector("zero", U().dimensions(), vector::zero));
    surfaceVectorField U_R (IOobject("U_r", mesh().time().timeName(), mesh()), mesh(), dimensionedVector("zero", U().dimensions(), vector::zero));
    U_L.replace(0, fvc::interpolate(U().component(0), pos_, "reconstruct(U)"));
    U_L.replace(1, fvc::interpolate(U().component(1), pos_, "reconstruct(U)"));
    U_L.replace(2, fvc::interpolate(U().component(2), pos_, "reconstruct(U)"));
    U_R.replace(0, fvc::interpolate(U().component(0), neg_, "reconstruct(U)"));
    U_R.replace(1, fvc::interpolate(U().component(1), neg_, "reconstruct(U)"));
    U_R.replace(2, fvc::interpolate(U().component(2), neg_, "reconstruct(U)"));

    tmp<surfaceScalarField> phi_L = U_L & mesh().Sf();
    tmp<surfaceScalarField> phi_R = U_R & mesh().Sf();

    // Flux relative to mesh movement
    if (mesh().moving())
    {
        fvc::makeRelative(phi_L.ref(), U());
        fvc::makeRelative(phi_R.ref(), U());
    }
    surfaceScalarField un_L = phi_L/mesh().magSf();
    surfaceScalarField un_R = phi_R/mesh().magSf();

    un_L -= MRFFaceVelocity();
    un_R -= MRFFaceVelocity();

    const volScalarField& T(thermo().T());
    tmp< volScalarField > E = thermo().he(p(), T)+0.5*magSqr(U());

    // Critical acoustic velocity (Liou 2006)
    tmp< volScalarField > gamma = thermo().gamma();
    tmp< volScalarField > H
    (
        (max(E(),dimensionedScalar("0", E().dimensions(), SMALL)) +
         max(p()/rho(),dimensionedScalar("0", p().dimensions()/rho().dimensions(), SMALL)))
    );
    H->rename("H");

    tmp< volScalarField > c = sqrt(2.0*(gamma()-1.0)/(gamma()+1.0)*H());
    c->rename("c");
    gamma.clear();

    tmp< surfaceScalarField > c_L = fvc::interpolate(c(), pos_, "reconstruct(T)");
	tmp< surfaceScalarField > c_R = fvc::interpolate(c(), neg_, "reconstruct(T)");

    c_L = sqr(c_L())/max(c_L(), un_L);
    c_R = sqr(c_R())/max(c_R(),-un_R);
    tmp< surfaceScalarField > c_face(min(c_L(),c_R()));
    c_L.clear();
    c_R.clear();

    // Critical Mach number
    tmp< surfaceScalarField > Mach_L(un_L/c_face());

    // Split Mach numbers
    tmp<surfaceScalarField> Mach_plus_L = 
        calcForEachFace
        (
            [](const scalar& MLf)
            {  
                if (mag(MLf) < 1.0)
                {
                    scalar ML2p =  0.25*sqr(MLf+1);
                    scalar ML2m = -0.25*sqr(MLf-1);
                    //return ML2p;             // beta = 0
                    return ML2p*(1 - 2*ML2m);  // beta = 1/8
                }
                else
                {
                    return max(MLf, 0);
                }
            },
            Mach_L()
        );

    // Pressure flux
    tmp<surfaceScalarField> p_plus_L =
        calcForEachFace
        (
            [](const scalar& MLf)
            {
                if (mag(MLf) < 1.0)
                {
                    scalar ML2p =  0.25*sqr(MLf+1);
                    scalar ML2m = -0.25*sqr(MLf-1);
                    return ML2p*(2 - MLf - 3*MLf*ML2m);  //alpha = 3/16
                }
                else
                {
                    return (MLf > 0 ? 1.0 : 0.0);
                }
            },
            Mach_L()
        );
    Mach_L.clear();

    tmp< surfaceScalarField > Mach_R(un_R/c_face());

    // Split Mach numbers
    tmp<surfaceScalarField> Mach_minus_R = 
        calcForEachFace
        (
            [](const scalar& MRf)
            {  
                if (mag(MRf) < 1.0)
                {
                    scalar MR2m = -0.25*sqr(MRf-1);
                    scalar MR2p =  0.25*sqr(MRf+1);
                    //return MR2m;             // beta = 0
                    return MR2m*(1 + 2*MR2p);  // beta = 1/8
                }
                else
                {
                    return min(MRf, 0);
                }
            },
            Mach_R()
        );

    // Pressure flux
    tmp<surfaceScalarField> p_minus_R =
        calcForEachFace
        (
            [](const scalar& MRf)
            {
                if (mag(MRf) < 1.0)
                {
                    scalar MR2m = -0.25*sqr(MRf-1);
                    scalar MR2p =  0.25*sqr(MRf+1);
                    return MR2m*(-2 - MRf + 3*MRf*MR2p);  //alpha = 3/16
                }
                else
                {
                    return (MRf < 0 ? 1.0 : 0.0);
                }
            },
            Mach_R()
        );
    Mach_R.clear();

    tmp< surfaceScalarField > p_L (fvc::interpolate(p(), pos_, "reconstruct(rho)"));
	tmp< surfaceScalarField > p_R (fvc::interpolate(p(), neg_, "reconstruct(rho)"));

    tmp< surfaceScalarField > rho_L (fvc::interpolate(rho(), pos_, "reconstruct(rho)"));
	tmp< surfaceScalarField > rho_R (fvc::interpolate(rho(), neg_, "reconstruct(rho)"));

    tmp< surfaceScalarField > Mach_1_2 = Mach_plus_L() + Mach_minus_R();
    Mach_plus_L.clear();
    Mach_minus_R.clear();

    surfaceScalarField p_1_2 = p_plus_L()*p_L() + p_minus_R()*p_R();

    // Low Mach number diffusive term
    tmp< surfaceScalarField > M_mean = 0.5*(sqr(un_L) + sqr(un_R))/sqr(c_face());   // Mean local Mach number
    tmp< surfaceScalarField > MDiff = -0.25*max((1.0-M_mean()),0.0)*(p_R() - p_L())/(0.5*(rho_L()+rho_R())*sqr(c_face()));    // 0 < K_p < 1, Liou suggest 0.25
    M_mean.clear();

    // Mach_1_2.ref() += MDiff();
    Mach_1_2.ref() +=
        calcForEachFace
        (
            [](const scalar& Mach_1_2, const scalar& MDiff)
            {
                if 
                (
                    (Mach_1_2 > 0.0 && Mach_1_2 + MDiff <= 0.0) ||
                    (Mach_1_2 < 0.0 && Mach_1_2 + MDiff >= 0.0)
                )
                {
                    return 0.2*MDiff;
                }
                else
                {
                    return MDiff;
                }
            },
            Mach_1_2(),
            MDiff()
        );

    MDiff.clear();

    // Low Mach number diffusive term
    if(lowMach_)
    {
        tmp< surfaceScalarField > pDiff = -0.25*p_plus_L()*p_minus_R()*(rho_L()+rho_R())*c_face()*(un_R-un_L); // 0 < Ku < 1; Liou suggests 0.75

        pDiff->setOriented(false);

        p_1_2 +=  pDiff();
        pDiff.clear();
    }

    Mach_1_2->setOriented(true);

    p_L.clear();
    p_R.clear();
    p_plus_L.clear();
    p_minus_R.clear();

    tmp<surfaceVectorField> U_f = surfaceFieldSelect(U_L, U_R, Mach_1_2(), 0);

    surfaceScalarField rhoa_LR  = Mach_1_2()*c_face()*surfaceFieldSelect(rho_L, rho_R, Mach_1_2(), 0);
    surfaceVectorField rhoaU_LR = rhoa_LR*U_f();

    tmp< surfaceScalarField > H_L (fvc::interpolate(H(), pos_, "reconstruct(T)"));
	tmp< surfaceScalarField > H_R (fvc::interpolate(H(), neg_, "reconstruct(T)"));

    surfaceScalarField rhoah_LR = rhoa_LR*(surfaceFieldSelect(H_L, H_R, Mach_1_2(), 0));

    U_f.clear();
    H.clear();

    phi = rhoa_LR*mesh().magSf();
    phiUp = rhoaU_LR*mesh().magSf() + p_1_2*mesh().Sf();
    phiEp = rhoah_LR*mesh().magSf() + p_1_2*MRFFaceVelocity()*mesh().magSf();

    if (mesh().moving())
    {
        phiEp += p_1_2 * fvc::meshPhi(U()); 
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
