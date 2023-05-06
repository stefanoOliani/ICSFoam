forAll(subTimeMeshes,K)
{
	// Reset pseudo Co number before every time step
	if (solnControl[0].dict().lookupOrDefault<bool>("resetPseudo", false) )
	{
		pseudoCoField[K] == solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
	}
}

scalar pseudoCoNumMin =
	solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNumMin", 0.1);

scalar pseudoCoNumMax =
	solnControl[0].dict().lookupOrDefault<scalar>("pseudoCoNumMax", 25);

scalar pseudoCoNumMaxIncr = solnControl[0].dict().lookupOrDefault<scalar>
(
	"pseudoCoNumMaxIncreaseFactor",
	1.25
);

scalar pseudoCoNumMinDecr = solnControl[0].dict().lookupOrDefault<scalar>
(
	"pseudoCoNumMinDecreaseFactor",
	0.1
);

HB.addingSource(contSource, scalarVarsGlob, 0);
HB.addingSource(momSource, vectorVarsGlob, 0, true);
HB.addingSource(energySource, scalarVarsGlob, 1);

scalar coNumRatio(0);

if (initRes.valid())
{
	if (!solnControl[0].firstIter() && prevRes.valid())
	{
		residualsIO& initialRes(initRes());
		residualsIO& previousRes(prevRes());

		scalar normInitSqr(0);
		scalar normPrevSqr(0);

		forAll(initialRes.vInitRes(), J)
		{
			normInitSqr += magSqr(initialRes.sInitRes()[2*J])
					     + magSqr(initialRes.sInitRes()[2*J+1])
						 + magSqr(initialRes.vInitRes()[J]);

			normPrevSqr += magSqr(previousRes.sInitRes()[2*J])
						 + magSqr(previousRes.sInitRes()[2*J+1])
						 + magSqr(previousRes.vInitRes()[J]);
		}

		scalar normInit = Foam::sqrt(normInitSqr);

		scalar normPrev = Foam::sqrt(normPrevSqr);

		coNumRatio = normPrev/normInit;

		coNumRatio = max(min(coNumRatio, pseudoCoNumMaxIncr), pseudoCoNumMinDecr);
	}

	prevRes.clear();
	prevRes.set(new residualsIO(initRes()));
}

forAll(subTimeMeshes,K)
{
	const fvMesh& mesh = subTimeMeshes[K];

	if (mesh.changing())
	{
		MRF[K].update();
		subTimeMeshes[K].setPhi().oldTime();
	}

	#include "setCoAndDeltaT.H"
}

PtrList<volScalarField> scalarVarsPrevIter(2*nO);
PtrList<volVectorField> vectorVarsPrevIter(nO);

forAll(vectorVarsGlob, i)
{
	scalarVarsPrevIter.set(2*i, new volScalarField(scalarVarsGlob[i][0].name() + "PrevIter", scalarVarsGlob[i][0]));
    vectorVarsPrevIter.set(i, new volVectorField(vectorVarsGlob[i][0].name() + "PrevIter", vectorVarsGlob[i][0]));
    scalarVarsPrevIter.set(2*i+1, new volScalarField(scalarVarsGlob[i][1].name() + "PrevIter", scalarVarsGlob[i][1]));
}

coupledMatrix eqSystem(subTimeMeshes[0], 2*nO, nO, true);

forAll(subTimeMeshes,K)
{
	const fvMesh& mesh = subTimeMeshes[K];

	volScalarField& rho_ = scalarVarsGlob[K][0];
	volVectorField& U_ = U[K];
	volVectorField& rhoU_ = vectorVarsGlob[K][0];
	volScalarField& rhoE_ = scalarVarsGlob[K][1];

	MRF[K].correctBoundaryVelocity(U_);

	convectiveFlux[K].MRFFaceVelocity().setOriented(true);
	convectiveFlux[K].MRFOmega() = MRF[K].omega()();

	surfaceVectorField phiUp
	(
		IOobject
		(
			"phiUp",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimArea*rhoU_.dimensions()*U_.dimensions()
	);

	surfaceScalarField phiEp
	(
		IOobject
		(
			"phiEp",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimArea*rhoE_.dimensions()*U_.dimensions()
	);

	convectiveFlux[K].calcFlux(phi[K], phiUp, phiEp);

	coupledMatrix eqSystemBlock(mesh, 2, 1, true); //true = delta form

	// Update residuals
	#include "residualsUpdate.H"

	scalarField ddtCoeff(mesh.V().size());

	// Set based on latest rPseudoDeltaT
	ddtCoeff = max(
					max(fvm::ddt(rho_)->diag(),fvm::ddt(rhoU_)->diag()),
					fvm::ddt(rhoE_)->diag()
				   )/mesh.V();

	convectiveFlux[K].createConvectiveJacobian(eqSystemBlock, ddtCoeff);

	if (!inviscid)
	{
		viscousFlux[K].createViscousJacobian(eqSystemBlock);
	}

	forAll(subTimeMeshes, J)
	{
		if (J==K)
		{
			HB.addBlock(eqSystemBlock.dSByS(0,0),J,K);
			HB.addBlock(eqSystemBlock.dVByV(0,0),J,K);
			HB.addBlock(eqSystemBlock.dSByS(1,1),J,K);

			eqSystem.dSByS(2*J, 2*K) += eqSystemBlock.dSByS(0,0);
			eqSystem.dVByV(J, K) += eqSystemBlock.dVByV(0,0);
			eqSystem.dSByS(2*J+1, 2*K+1) += eqSystemBlock.dSByS(1,1);

			eqSystem.dSByS(2*J, 2*K+1) = eqSystemBlock.dSByS(0,1);
			eqSystem.dSByS(2*J+1, 2*K) = eqSystemBlock.dSByS(1,0);

			eqSystem.dSByV(2*J, K) = eqSystemBlock.dSByV(0,0);
			eqSystem.dSByV(2*J+1, K) = eqSystemBlock.dSByV(1,0);

			eqSystem.dVByS(J, 2*K) = eqSystemBlock.dVByS(0,0);
			eqSystem.dVByS(J, 2*K+1) = eqSystemBlock.dVByS(0,1);
		}
		else
		{
			blockFvMatrix<scalar,scalar> HBRho(subTimeMeshes[J]);
			blockFvMatrix<vector,tensor> HBRhoU(subTimeMeshes[J]);
			blockFvMatrix<scalar,scalar> HBRhoE(subTimeMeshes[J]);
			HB.addBlock(HBRho,J,K);
			HB.addBlock(HBRhoU,J,K);
			HB.addBlock(HBRhoE,J,K);

			//To avoid unallocated diagonal pointer in case
			// HB zones are not active
			eqSystem.dSByS(2*J, 2*K).diag() = Zero;
			eqSystem.dVByV(J, K).diag() = Zero;
			eqSystem.dSByS(2*J+1, 2*K+1).diag() = Zero;

			eqSystem.dSByS(2*J, 2*K) = HBRho;
			eqSystem.dVByV(J, K) = HBRhoU;
			eqSystem.dSByS(2*J+1, 2*K+1) = HBRhoE;

			eqSystem.dSByS(2*J, 2*K+1).diag() = Zero;
			eqSystem.dSByS(2*J+1, 2*K).diag() = Zero;

			eqSystem.dSByV(2*J, K).diag() = Zero;
			eqSystem.dSByV(2*J+1, K).diag() = Zero;

			eqSystem.dVByS(J, 2*K).diag() = Zero;
			eqSystem.dVByS(J, 2*K+1).diag() = Zero;
		}
	}
}

#include "setUpSolIncrements.H"

initRes.clear();

initRes.set(new residualsIO(eqSystem.solveForIncr
		(
			scalarVarsPrevIter,
			vectorVarsPrevIter,
			scalarSolIncrement,
			vectorSolIncrement
		)));

forAll(subTimeMeshes,K)
{
	#include "setUpResiduals.H"

	#include "boundLocalTimeStep.H"

	turbulence[K].correct();

	#include "updateFields.H"
}
