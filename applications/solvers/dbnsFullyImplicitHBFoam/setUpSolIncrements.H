PtrList<volScalarField> scalarSolIncrement(2*nO);
PtrList<volVectorField> vectorSolIncrement(nO);

forAll(vectorSolIncrement, i)
{
	scalarSolIncrement.set
	(
		2*i,
		new volScalarField
		(
			IOobject
			(
				"rhoIncr",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i],
			dimensionedScalar("rhoIncr", dimDensity, 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);

	vectorSolIncrement.set
	(
		i,
		new volVectorField
		(
			IOobject
			(
				"rhoUIncr",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i],
			dimensionedVector("rhoUIncr", dimDensity*dimVelocity, Zero),
			zeroGradientFvPatchScalarField::typeName
		)
	);

	scalarSolIncrement.set
	(
		2*i+1,
		new volScalarField
		(
			IOobject
			(
				"rhoEIncr",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i],
			dimensionedScalar("rhoEIncr", dimDensity*dimVelocity*dimVelocity, 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}
