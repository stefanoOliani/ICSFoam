Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermoMesh
(
    psiThermo::New(mesh)
);

psiThermo& thermoMesh = pThermoMesh();

volScalarField& p = thermoMesh.p();
volScalarField& e = thermoMesh.he();
volScalarField& T = thermoMesh.T();

if (e.name() != "e")
{
    FatalErrorInFunction
        << "Only energy type internalEnergy supported."
        << nl << exit(FatalError);
}

const volScalarField& mu = thermoMesh.mu();

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

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermoMesh.rho()
);

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
	mesh.Sf() & fvc::interpolate(thermoMesh.rho()*U)
);



PtrList<psiThermo> pThermoHB(nO);
PtrList<compressible::turbulenceModel> turbulenceHB(nO);

PtrList<surfaceScalarField> phiHB(nO);
PtrList<volVectorField> UHB(nO);
PtrList<volScalarField> pHB(nO);
PtrList<volScalarField> THB(nO);
PtrList<volScalarField> rhoHB(nO);


Info<< "Reading thermophysical properties\n" << endl;

forAll(pThermoHB,i)
{
	autoPtr<psiThermo> PThermo
	(
	    psiThermo::New(subTimeMeshes[i])
	);

	pThermoHB.set
	(
		i,
		PThermo
	);
}


Info<< "Reading field UHB" << nl << endl;

forAll(UHB,i)
{
	UHB.set
	(
		i,
	    new volVectorField
	    (
	        IOobject
	        (
	            "U",
	            runTime.timeName(),
				subTimeMeshes[i],
	            IOobject::MUST_READ,
	            IOobject::NO_WRITE
	        ),
			subTimeMeshes[i]
	    )
	);

	pHB.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"p",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i]
		)
	);

	THB.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"T",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i]
		)
	);

	rhoHB.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"rho",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i]
		)
	);
}


forAll(phiHB,i)
{
	volVectorField rhoUHB = rhoHB[i]*UHB[i];

	phiHB.set
	(
		i,
	    new surfaceScalarField
	    (
	        IOobject
	        (
	            "phi",
	            runTime.timeName(),
				subTimeMeshes[i],
	            IOobject::READ_IF_PRESENT,
	            IOobject::NO_WRITE
	        ),
			subTimeMeshes[i].Sf() & fvc::interpolate(rhoUHB)
	    )
	);
}

Info<< "Creating turbulence HB model\n" << endl;
forAll(turbulenceHB,i)
{
	autoPtr<compressible::turbulenceModel> Turbulence
	(
			compressible::turbulenceModel::New
			(
				rhoHB[i],
				UHB[i],
				phiHB[i],
				pThermoHB[i]
			)
	);

	turbulenceHB.set(i, Turbulence);
}

