// Construct cont eqn: [A] \Delta W = [R]
    volScalarField rhoR =
(
    -fvc::div(phi[K])
);

// Construct mom eqn: [A] \Delta W = [R]
 volVectorField rhoUR =
(
    -fvc::div(phiUp)
);

// Construct energ eqn: [A] \Delta W = [R]
volScalarField rhoER =
(
    -fvc::div(phiEp)
);

// Viscous flow
if (!inviscid)
{
    volScalarField muEff("muEff", turbulence[K].muEff());

    volScalarField alphaEff("alphaEff", turbulence[K].alphaEff());

    volTensorField tauMC = muEff*dev2(Foam::T(fvc::grad(U_)));

    rhoUR += fvc::laplacian(muEff,U_);
    rhoUR += fvc::div(tauMC);

    surfaceVectorField sigmaDotU
    (
        "sigmaDotU",
        (
          fvc::interpolate(muEff)*fvc::interpolate(fvc::grad(U_))
          + fvc::interpolate(tauMC)
        )
        & fvc::interpolate(U_)
    );

    rhoER += fvc::div(sigmaDotU & mesh.Sf());

    volScalarField eCalc("eCalc", rhoE_/rho_-0.5*magSqr(U_));

    rhoER += fvc::laplacian(alphaEff, eCalc, "laplacian(alphaEff,e)");
}

// Pseudo and real-time contribution

tmp<fvScalarMatrix> ddtRho = fvm::ddt(rho_);
tmp<fvVectorMatrix> ddtRhoU = fvm::ddt(rhoU_);
tmp<fvScalarMatrix> ddtRhoE = fvm::ddt(rhoE_);

// Work around bug causing timeIndex not to be stored for ddt0 field (CrankNicolson)
// Seems to be caused by commented out line in GeometricField::storeOldTimes(),
if (mesh.foundObject<volScalarField>("ddt0(rho)"))
{
    const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>("ddt0(rho)")
    ).timeIndex() = mesh.time().timeIndex();
}

if (mesh.foundObject<volVectorField>("ddt0(rhoU)"))
{
    const_cast<volVectorField&>
    (
        mesh.lookupObject<volVectorField>("ddt0(rhoU)")
    ).timeIndex() = mesh.time().timeIndex();
}

if (mesh.foundObject<volScalarField>("ddt0(rhoE)"))
{
    const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>("ddt0(rhoE)")
    ).timeIndex() = mesh.time().timeIndex();
}

eqSystemBlock.dSByS(0,0).source() = rhoR.primitiveField()*mesh.V() + contSource[K];
eqSystemBlock.dVByV(0,0).source() = rhoUR.primitiveField()*mesh.V() + momSource[K];
eqSystemBlock.dSByS(1,1).source() = rhoER.primitiveField()*mesh.V() + energySource[K];
