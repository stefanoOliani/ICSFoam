PtrList<volScalarField> scalarVars(2);
PtrList<volVectorField> vectorVars(1);

PtrList<volScalarField> scalarVarsScalarTransp(1);

bool scalarTransport = false;
dictionary scalarDict("fvSolution");

if (mesh.solutionDict().findDict("scalarTransport"))
{
	scalarDict = mesh.solutionDict().subDict("scalarTransport");
	scalarTransport = scalarDict.getOrDefault("active", false);
}

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);

psiThermo& thermo = pThermo();

volScalarField& p = thermo.p();
volScalarField& e = thermo.he();
volScalarField& T = thermo.T();

if (e.name() != "e")
{
    FatalErrorInFunction
        << "Only energy type internalEnergy supported."
        << nl << exit(FatalError);
}

const volScalarField& mu = thermo.mu();

bool inviscid(false);

if (gMax(mu.internalField()) > 0.0)
{
    Info << "Viscous analysis detected" << nl << endl;
}
else
{
    inviscid = true;
    Info << "Inviscid analysis detected" << nl << endl;
}

Info<< "Reading field U" << nl << endl;

volVectorField U
(
	IOobject
	(
		"U",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volVectorField URel
(
	IOobject
	(
		"URel",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	U
);

scalarVars.set
(
    0,
    new volScalarField
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    )
);

volScalarField& rho = scalarVars[0];

vectorVars.set
(
    0,
    new volVectorField
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        rho*U
    )
);

volVectorField& rhoU = vectorVars[0];

scalarVars.set
(
    1,
    new volScalarField
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        rho*(e + 0.5*magSqr(U))
    )
);

volScalarField& rhoE = scalarVars[1];

surfaceScalarField phi
(
	IOobject
	(
		"phi",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh.Sf() & fvc::interpolate(rhoU)
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating flux scheme\n" << endl;
autoPtr<convectiveFluxScheme> pFlux
(
	convectiveFluxScheme::New
	(
		mesh.schemesDict(),
		thermo,
		rho,
		U,
		p
	)
);

convectiveFluxScheme& flux = pFlux();

Info<< "Creating viscous flux scheme\n" << endl;
viscousFluxScheme viscousFlux(mesh.schemesDict(), thermo, rho, U, p, turbulence);

volScalarField rPseudoDeltaT
(
	IOobject
	(
		"rPseudoDeltaT",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("rPseudoDeltaT", dimless/dimTime, 0.0),
	zeroGradientFvPatchScalarField::typeName
);

scalarField ddtCoeff(mesh.V().size());

autoPtr<volTensorField> pTauMC(nullptr);

if (!inviscid)
{
    // tauMC is stored for the benefit of the maxwellSlip BC - maybe better
    // to modify the BC to recalculate it rather
    pTauMC.set
    (
        new volTensorField
        (
            IOobject
            (
                "tauMC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            turbulence->muEff()*dev2(Foam::T(fvc::grad(U)))
        )
    );
}


volScalarField pseudoCoField
(
	IOobject
	(
		"pseudoCoField",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar
	(
		"pseudoCoNumInit",
		dimless,
		solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
	)
);

uniformDimensionedScalarField pseudoCoNum
(
	IOobject
	(
		"pseudoCoNum",
		runTime.timeName(),
		"uniform",
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	dimensionedScalar
	(
		"pseudoCoNumInit",
		dimless,
		solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
	)
 );

//Store previous residuals to control pseudoCo increase
autoPtr<residualsIO> prevRes(nullptr);
autoPtr<residualsIO> initRes(nullptr);

autoPtr<residualsIO> initResScalarTransp(nullptr);

IOMRFCoupledZoneList MRF(mesh);
MRF.correctBoundaryVelocity(U);

IOMRFTranslatingZoneList MRFTrans(mesh);
MRFTrans.correctBoundaryVelocity(U);

MRF.makeRelative(URel);
MRFTrans.makeRelative(URel);

autoPtr<volScalarField> sPtr(nullptr);

if (scalarTransport)
{
	Info << "Scalar transport equation detected" << nl << endl;
	Info << "Reading s field" << nl << endl;

	sPtr.set
	(
		new volScalarField
		(
			IOobject
			(
				"s",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		)
	);

	const volScalarField& s = sPtr();

	scalarVarsScalarTransp.set
	(
		0,
		new volScalarField
		(
			IOobject
			(
				"rhoS",
				runTime.timeName(),
				mesh,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE
			),
			rho*s
		)
	);
}

#include "createFvOptions.H"
