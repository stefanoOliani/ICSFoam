PtrList<volScalarField> scalarSolIncrement(2);
PtrList<volVectorField> vectorSolIncrement(1);

scalarSolIncrement.set
(
	0,
	new volScalarField
	(
		IOobject
		(
			"rhoIncr",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("rhoIncr", dimDensity, 0.0),
		zeroGradientFvPatchScalarField::typeName
	)
);

vectorSolIncrement.set
(
	0,
	new volVectorField
	(
		IOobject
		(
			"rhoUIncr",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedVector("rhoUIncr", dimDensity*dimVelocity, Zero),
		zeroGradientFvPatchScalarField::typeName
	)
);

scalarSolIncrement.set
(
	1,
	new volScalarField
	(
		IOobject
		(
			"rhoEIncr",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("rhoEIncr", dimDensity*dimVelocity*dimVelocity, 0.0),
		zeroGradientFvPatchScalarField::typeName
	)
);

