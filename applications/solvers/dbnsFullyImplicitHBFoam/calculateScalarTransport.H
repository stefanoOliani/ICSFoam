HB.addingSource(scalarTranspSource, scalarVarsGlobScalarTransp, 0);

PtrList<volScalarField> scalarVarsPrevIterScalarTransp(nO);
PtrList<volVectorField> dummyVectorVarsPrevIter(0);

forAll(scalarVarsGlobScalarTransp, i)
{
	scalarVarsPrevIterScalarTransp.set
	(
		i,
		new volScalarField(scalarVarsGlobScalarTransp[i][0].name() + "PrevIter", scalarVarsGlobScalarTransp[i][0])
	);
}

coupledMatrix scalarTranspEqSystem(subTimeMeshes[0], nO, 0, true);

forAll(subTimeMeshes,K)
{
	const fvMesh& mesh = subTimeMeshes[K];

	volScalarField& rhoS_ = scalarVarsGlobScalarTransp[K][0];
	volScalarField& s_ = s[K];
	volScalarField& rho_ = scalarVarsGlob[K][0];
	volVectorField& U_ = U[K];

	surfaceScalarField& MRFFaceVelocity = convectiveFlux[K].MRFFaceVelocity();

	surfaceScalarField phiSp
	(
		IOobject
		(
			"phiSp",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimArea*rhoS_.dimensions()*U_.dimensions()
	);

	surfaceScalarField pos_
	(
	    IOobject
	    (
	        "pos",
	        runTime.timeName(),
	        mesh
	    ),
	    mesh,
	    dimensionedScalar("pos", dimless, 1.0)
	);

	surfaceScalarField neg_
	(
	    IOobject
	    (
	        "neg",
	        runTime.timeName(),
	        mesh
	    ),
	    mesh,
	    dimensionedScalar("neg", dimless, -1.0)
	);

	surfaceVectorField n = mesh.Sf()/mesh.magSf();

	tmp< surfaceScalarField > s_l = fvc::interpolate(s_, pos_, "reconstruct(s)");
	tmp< surfaceScalarField > s_r = fvc::interpolate(s_, neg_, "reconstruct(s)");

	tmp< surfaceScalarField > rho_l (fvc::interpolate(rho_, pos_, "reconstruct(rho)"));
	tmp< surfaceScalarField > rho_r (fvc::interpolate(rho_, neg_, "reconstruct(rho)"));

	//Conservative variables reconstruction at interface
	surfaceScalarField rhoS_l(rho_l() * s_l());
	surfaceScalarField rhoS_r(rho_r() * s_r());

	surfaceVectorField U_l (IOobject("U_l", mesh.time().timeName(), mesh), mesh, dimensionedVector("zero", U_.dimensions(), vector::zero));
	surfaceVectorField U_r (IOobject("U_r", mesh.time().timeName(), mesh), mesh, dimensionedVector("zero", U_.dimensions(), vector::zero));
	U_l.replace(0, fvc::interpolate(U_.component(0), pos_, "reconstruct(U)"));
	U_l.replace(1, fvc::interpolate(U_.component(1), pos_, "reconstruct(U)"));
	U_l.replace(2, fvc::interpolate(U_.component(2), pos_, "reconstruct(U)"));
	U_r.replace(0, fvc::interpolate(U_.component(0), neg_, "reconstruct(U)"));
	U_r.replace(1, fvc::interpolate(U_.component(1), neg_, "reconstruct(U)"));
	U_r.replace(2, fvc::interpolate(U_.component(2), neg_, "reconstruct(U)"));

	tmp< surfaceScalarField > magU_l = U_l & n;
	tmp< surfaceScalarField > magU_r = U_r & n;

	// Roe averages
	dimensionedScalar rho0("rho0",dimDensity,VSMALL);
	tmp< surfaceScalarField > coefR = sqrt(max(rho0,rho_r())/max(rho0,rho_l()));  // CHECK: Is this form of coefR correct?

	surfaceVectorField RoeVelocity = (coefR()*U_r + U_l)/(coefR() + 1.0);
	surfaceScalarField  uProjRoe = RoeVelocity & n;

	if (mesh.moving())
	{
		surfaceScalarField meshVelocity = fvc::meshPhi(U_)()/mesh.magSf();

		uProjRoe -= meshVelocity;
	}
	uProjRoe -= MRFFaceVelocity;

	surfaceScalarField Lambda = mag(uProjRoe);

	//Harten's entropy correction
	surfaceScalarField epsilon =  0.05* Lambda;

	forAll(epsilon, facei)
	{
	    if(Lambda[facei] < epsilon[facei])
	    {
	        Lambda[facei] = (sqr(Lambda[facei]) + sqr(epsilon[facei]))/(2.0*epsilon[facei]);
	    }
	}

	forAll(epsilon.boundaryField(), patchi)
	{
	    const scalarField& pEpsilon = epsilon.boundaryField()[patchi];

	    scalarField& pLambda = Lambda.boundaryFieldRef()[patchi];

	    forAll(pEpsilon, facei)
	    {
	    	 if(pLambda[facei] < pEpsilon[facei])
	    	 {
	    		 pLambda[facei] = (sqr(pLambda[facei]) + sqr(pEpsilon[facei]))/(2.0*pEpsilon[facei]);
	    	 }
	    }
	}

	// Differences of conservative variables at interface
	surfaceScalarField diffRhoS = rhoS_r - rhoS_l;

	phiSp = -0.5*mesh.magSf()*Lambda*diffRhoS;
	phiSp.setOriented(true);

	surfaceScalarField rhoUNorm_l = rho_l()*magU_l();
	surfaceScalarField rhoUNorm_r = rho_r()*magU_r();

	phiSp += 0.5*mesh.magSf()*(rhoUNorm_l*s_l() + rhoUNorm_r*s_r());

	if (mesh.moving())
	{
		surfaceScalarField meshVelocity = fvc::meshPhi(U_)/mesh.magSf();

		phiSp -= 0.5*mesh.magSf()*meshVelocity*(rhoS_l + rhoS_r);
	}

	phiSp -= 0.5*mesh.magSf()*MRFFaceVelocity*(rhoS_l + rhoS_r);

	s_l.clear();
	s_r.clear();

	coupledMatrix scalarTranspEqSystemBlock(mesh, 1, 0, true); //true = delta form

	// Update residuals
    // Construct scalar transport eqn: [A] \Delta W = [R]
	volScalarField rhoSR =
    (
        -fvc::div(phiSp)
    );

    // Viscous flow
    if (!inviscid)
    {
    	scalar D_(0.0);
    	D_ = scalarDict.get<scalar>("D");

        volScalarField D
        (
            IOobject
            (
                "D",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("D", phi[K].dimensions()/dimLength, D_)
        );

        rhoSR += fvc::laplacian(D,s_);
    }

    // Pseudo and real-time contribution
    tmp<fvScalarMatrix> ddtRhoS = fvm::ddt(rhoS_);

    if (mesh.foundObject<volScalarField>("ddt0(rhoS)"))
    {
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>("ddt0(rhoS)")
        ).timeIndex() = mesh.time().timeIndex();
    }

    scalarTranspEqSystemBlock.dSByS(0,0).source() = rhoSR.primitiveField()*mesh.V() + scalarTranspSource[K];

    surfaceScalarField w(mesh.weights());
    w.setOriented(false);

    surfaceScalarField projU_l = U_l & n;
    surfaceScalarField projU_r = U_r & n;

	blockFvMatrix<scalar,scalar>& dScalarFluxByRhoS = scalarTranspEqSystemBlock.dSByS(0,0);

	dScalarFluxByRhoS.insertBlock(projU_l,projU_r);

	// Moving mesh part
	if (mesh.moving())
	{
		const surfaceScalarField& meshPhi(mesh.phi());
		tmp<blockFvMatrix<scalar,scalar> > divMeshPhi = fvj::div(w, meshPhi);

		dScalarFluxByRhoS -= divMeshPhi();
	}

	tmp<blockFvMatrix<scalar,scalar> > MRFdivMeshPhi = fvj::div(w, MRFFaceVelocity*mesh.magSf());

	dScalarFluxByRhoS -= MRFdivMeshPhi();

	MRFdivMeshPhi.clear();

	surfaceScalarField lambdaConv
	(
		IOobject
		(
			"lambdaConv",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("lambdaConv", dimVelocity, 0.0)
	);

	if (mesh.moving())
	{
		lambdaConv = (mag
				(
					(fvc::interpolate(U_) & mesh.Sf()/mesh.magSf())
					- MRFFaceVelocity
					- fvc::meshPhi(U_)/mesh.magSf()
				));
	}
	else
	{
		lambdaConv = mag((fvc::interpolate(U_) & mesh.Sf()/mesh.magSf()) - MRFFaceVelocity);
	}

	dScalarFluxByRhoS.insertDissipationBlock(lambdaConv);

    // Calculate diagonal contribution of the geometric boundaries

    forAll(mesh.boundaryMesh(), patchi)
    {
        const fvPatch& patch = mesh.boundary()[patchi];

        if(!patch.coupled() && !isA<emptyFvPatch>(patch))
        {
            tmp<scalarField> dScalarTranspFluxds;

            const volScalarField::Boundary& sbf = s_.boundaryField();
            const volVectorField::Boundary& ubf = U_.boundaryField();

            const scalarField& wP = mesh.weights().boundaryField()[patchi];

            // ValueInternalCoeffs
            // fixedValue   -> 0
            // zeroGradient -> 1

            tmp < scalarField > stVIC = sbf[patchi].valueInternalCoeffs(wP);
            const scalarField& sVIC = stVIC();

            const vectorField& SfB = mesh.Sf().boundaryField()[patchi];

            const scalarField& rhoB = rho_.boundaryField()[patchi];
            tmp<vectorField> UB = ubf[patchi];

            tmp<scalarField> UrelBdotSf = UB() & SfB;

            if (mesh.moving())
            {
                const surfaceScalarField& meshPhi(mesh.phi());
                UrelBdotSf.ref() -= meshPhi.boundaryField()[patchi];
            }

            UrelBdotSf.ref() -= MRFFaceVelocity.boundaryField()[patchi]*mesh.magSf().boundaryField()[patchi];

            dScalarTranspFluxds = rhoB * UrelBdotSf() * sVIC;

            tmp<scalarField> rhoI = rho_.boundaryField()[patchi].patchInternalField();
            tmp<scalarField> dsdRhoS = 1.0/rhoI();
            scalarField& diagScalarFluxByRho = dScalarFluxByRhoS.diag();

            forAll(patch, bfacei)
            {
                label iIntCell = patch.faceCells()[bfacei];

                diagScalarFluxByRho[iIntCell] += dScalarTranspFluxds()[bfacei] * dsdRhoS()[bfacei];
            }
        }
    }

    // Add viscous terms
	if (!inviscid)
	{
    	scalar D_(0.0);
    	D_ = scalarDict.get<scalar>("D");

        volScalarField D
        (
            IOobject
            (
                "D",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("D", phi[K].dimensions()/dimLength, D_)
        );

		surfaceScalarField lambdaVisc = fvc::interpolate(D)/fvc::interpolate(rho_);

		tmp<blockFvMatrix<scalar,scalar> > stab = fvj::laplacian(0.5*lambdaVisc, geometricOneField());

		dScalarFluxByRhoS -= stab();
	}

    // Add temporal terms
    // Temporal contribution (Weak form)
	scalarField ddtCoeff(mesh.V().size());

	// Set based on latest rPseudoDeltaT
	ddtCoeff = (fvm::ddt(rhoS_)->diag())/mesh.V();

    scalarField diagCoeff = ddtCoeff*mesh.V();
    dScalarFluxByRhoS.diag() += diagCoeff;

	forAll(subTimeMeshes, J)
	{
		if (J==K)
		{
			HB.addBlock(scalarTranspEqSystemBlock.dSByS(0,0),J,K);

			scalarTranspEqSystem.dSByS(J, K) += scalarTranspEqSystemBlock.dSByS(0,0);
		}
		else
		{
			blockFvMatrix<scalar,scalar> HBRhoS(subTimeMeshes[J]);

			HB.addBlock(HBRhoS,J,K);

			//To avoid unallocated diagonal pointer in case
			// HB zones are not active
			scalarTranspEqSystem.dSByS(J, K).diag() = Zero;

			scalarTranspEqSystem.dSByS(J, K) = HBRhoS;
		}
	}
}

PtrList<volScalarField> scalarSolIncrementScalarTransp(nO);
PtrList<volVectorField> dummyVectorSolIncrement(0);

forAll(scalarSolIncrementScalarTransp, i)
{
	scalarSolIncrementScalarTransp.set
	(
		i,
		new volScalarField
		(
			IOobject
			(
				"rhoSIncr",
				runTime.timeName(),
				subTimeMeshes[i],
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			subTimeMeshes[i],
			dimensionedScalar("rhoSIncr", dimDensity, 0.0),
			zeroGradientFvPatchScalarField::typeName
		)
	);
}

initResScalarTransp.clear();

initResScalarTransp.set(new residualsIO(scalarTranspEqSystem.solveForIncr
		(
			scalarVarsPrevIterScalarTransp,
			dummyVectorVarsPrevIter,
			scalarSolIncrementScalarTransp,
			dummyVectorSolIncrement
		)));

forAll(subTimeMeshes,K)
{
	volScalarField& rhoS_ = scalarVarsGlobScalarTransp[K][0];

	volScalarField& s_ = s[K];
	volScalarField& rho_ = scalarVarsGlob[K][0];
	volScalarField& rhoSIncr = scalarSolIncrementScalarTransp[K];
	rhoS_.ref() += rhoSIncr.internalField();

	// Compute internal fields
	s_.ref() = rhoS_.internalField()/rho_.internalField();

	// Correct boundary fields
	s_.correctBoundaryConditions();
	rhoS_.boundaryFieldRef() = rho_.boundaryField()*s_.boundaryField();
}
