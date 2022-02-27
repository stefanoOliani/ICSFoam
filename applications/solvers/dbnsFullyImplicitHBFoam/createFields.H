PtrList<PtrList<volScalarField>> scalarVarsGlob(nO);
PtrList<PtrList<volVectorField>> vectorVarsGlob(nO);

PtrList<scalarField> contSource(nO);
PtrList<vectorField> momSource(nO);
PtrList<scalarField> energySource(nO);

PtrList<psiThermo> pThermo(nO);
PtrList<compressible::turbulenceModel> turbulence(nO);
PtrList<convectiveFluxScheme> convectiveFlux(nO);
PtrList<viscousFluxScheme> viscousFlux(nO);

PtrList<surfaceScalarField> phi(nO);

PtrList<volVectorField> U(nO);

PtrList<volScalarField> pseudoCoField(nO);
PtrList<volScalarField> rPseudoDeltaT(nO);

autoPtr<residualsIO> prevRes(nullptr);
autoPtr<residualsIO> initRes(nullptr);

PtrList<IOMRFCoupledZoneList> MRF(nO);

PtrList<PtrList<volScalarField>> scalarVarsGlobScalarTransp(nO);
PtrList<volScalarField> s(nO);
PtrList<scalarField> scalarTranspSource(nO);

autoPtr<residualsIO> initResScalarTransp(nullptr);

bool scalarTransport = false;
dictionary scalarDict("fvSolution");

if (subTimeMeshes[0].solutionDict().findDict("scalarTransport"))
{
	scalarDict = subTimeMeshes[0].solutionDict().subDict("scalarTransport");
	scalarTransport = scalarDict.getOrDefault("active", false);
}

Info<< "Reading field U" << nl << endl;

forAll(U,i)
{
	U.set
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
	            IOobject::AUTO_WRITE
	        ),
			subTimeMeshes[i]
	    )
	);
}

Info<< "Reading thermophysical properties\n" << endl;

bool inviscid(true);

forAll(pThermo,i)
{
	autoPtr<psiThermo> PThermo
	(
	    psiThermo::New(subTimeMeshes[i])
	);

	psiThermo& thermo = PThermo();

	volScalarField& e = thermo.he();

	pThermo.set
	(
		i,
		PThermo
	);

	if (e.name() != "e")
	{
	    FatalErrorInFunction
	        << "Only energy type internalEnergy supported."
	        << nl << exit(FatalError);
	}

	const volScalarField& mu = thermo.mu();

	if (gMax(mu.internalField()) > 0.0)
	{
		inviscid = false;
	}
}

if (!inviscid)
{
    Info << "Viscous analysis detected" << nl << endl;
}
else
{
    Info << "Inviscid analysis detected" << nl << endl;
}

forAll(scalarVarsGlob, i)
{
	psiThermo& thermo = pThermo[i];

	scalarVarsGlob.set
	(
	    i,
	    new PtrList<volScalarField>(2)
	);

	vectorVarsGlob.set
	(
	    i,
	    new PtrList<volVectorField>(1)
	);

	scalarVarsGlob[i].set
	(
		0,
		new volScalarField
		(
			IOobject
			(
				"rho",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			thermo.rho()
		)
	);

	volScalarField& rho = scalarVarsGlob[i][0];
	volVectorField& Ui = U[i];

	vectorVarsGlob[i].set
	(
	    0,
	    new volVectorField
	    (
	        IOobject
	        (
	            "rhoU",
	            runTime.timeName(),
				subTimeMeshes[i],
	            IOobject::READ_IF_PRESENT,
	            IOobject::NO_WRITE
	        ),
	        rho*Ui
	    )
	);

	volScalarField& e = thermo.he();

	scalarVarsGlob[i].set
	(
	    1,
	    new volScalarField
	    (
	        IOobject
	        (
	            "rhoE",
	            runTime.timeName(),
				subTimeMeshes[i],
	            IOobject::READ_IF_PRESENT,
	            IOobject::NO_WRITE
	        ),
	        rho*(e + 0.5*magSqr(Ui))
	    )
	);

}

forAll(phi,i)
{
	volVectorField& rhoU = vectorVarsGlob[i][0];

	phi.set
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
	            IOobject::AUTO_WRITE
	        ),
			subTimeMeshes[i].Sf() & fvc::interpolate(rhoU)
	    )
	);
}

Info<< "Creating turbulence model\n" << endl;
forAll(turbulence,i)
{
	volScalarField& rho = scalarVarsGlob[i][0];

	autoPtr<compressible::turbulenceModel> Turbulence
	(
			compressible::turbulenceModel::New
			(
				rho,
				U[i],
				phi[i],
				pThermo[i]
			)
	);

	turbulence.set(i, Turbulence);
}

Info<< "Creating flux scheme\n" << endl;
forAll(convectiveFlux,i)
{
	volScalarField& rho = scalarVarsGlob[i][0];
	volVectorField& Ui = U[i];
	volScalarField& p = pThermo[i].p();

	autoPtr<convectiveFluxScheme> CFlux
	(
			convectiveFluxScheme::New
			(
				subTimeMeshes[0].schemesDict(),
				pThermo[i],
				rho,
				Ui,
				p
			)
	);

	convectiveFlux.set(i, CFlux);
}


Info<< "Creating viscous flux scheme\n" << endl;
forAll(viscousFlux,i)
{
	volScalarField& rho = scalarVarsGlob[i][0];
	volVectorField& Ui = U[i];
	volScalarField& p = pThermo[i].p();

	viscousFlux.set
	(
		i,
	    new viscousFluxScheme
	    (
			subTimeMeshes[0].schemesDict(),
			pThermo[i],
			rho,
			Ui,
			p,
			turbulence[i]
	    )
	);
}

forAll(rPseudoDeltaT,i)
{
	rPseudoDeltaT.set
	(
		i,
	    new volScalarField
		(
			IOobject
			(
				"rPseudoDeltaT",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i],
			dimensionedScalar("rPseudoDeltaT", dimless/dimTime, 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

forAll(pseudoCoField,i)
{
	pseudoCoField.set
	(
		i,
	    new volScalarField
		(
			IOobject
			(
				"pseudoCoField",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::READ_IF_PRESENT,
				IOobject::AUTO_WRITE
			),
			subTimeMeshes[i],
			dimensionedScalar
			(
				"pseudoCoNumInit",
				dimless,
				solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNum", 1)
			)
		)
	);
}

forAll(MRF,i)
{
	MRF.set
	(
		i,
	    new IOMRFCoupledZoneList(subTimeMeshes[i])
	);

	MRF[i].correctBoundaryVelocity(U[i]);
}


forAll(contSource,i)
{
	contSource.set
	(
		i,
	    new scalarField(U[i].primitiveField().size(), Zero)
	);

	momSource.set
	(
		i,
	    new vectorField(U[i].primitiveField().size(), Zero)
	);

	energySource.set
	(
		i,
	    new scalarField(U[i].primitiveField().size(), Zero)
	);
}


if (scalarTransport)
{
	Info << "Scalar transport equation detected" << nl << endl;
	Info << "Reading s field" << nl << endl;

	forAll(scalarVarsGlobScalarTransp, i)
	{
		s.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"s",
					runTime.timeName(),
					subTimeMeshes[i],
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				subTimeMeshes[i]
			)
		);

		scalarVarsGlobScalarTransp.set
		(
			i,
			new PtrList<volScalarField>(1)
		);

		volScalarField& rho = scalarVarsGlob[i][0];
		volScalarField& si = s[i];

		scalarVarsGlobScalarTransp[i].set
		(
		    0,
		    new volScalarField
		    (
		        IOobject
		        (
		            "rhoS",
		            runTime.timeName(),
					subTimeMeshes[i],
		            IOobject::READ_IF_PRESENT,
		            IOobject::NO_WRITE
		        ),
		        rho*si
		    )
		);
	}

	forAll(scalarTranspSource,i)
	{
		scalarTranspSource.set
		(
			i,
			new scalarField(U[i].primitiveField().size(), Zero)
		);
	}
}
