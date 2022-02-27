// Switched Evolution Relaxation
if (initRes.valid())
{
	if (!solnControl.firstIter() && prevRes.valid())
	{
		residualsIO& initialRes(initRes());
		residualsIO& previousRes(prevRes());

		scalar normInit = Foam::sqrt(magSqr(initialRes.sInitRes()[0]) +
							magSqr(initialRes.sInitRes()[1]) +
							magSqr(initialRes.vInitRes()[0]));

		scalar normPrev = Foam::sqrt(magSqr(previousRes.sInitRes()[0]) +
							magSqr(previousRes.sInitRes()[1]) +
							magSqr(previousRes.vInitRes()[0]));

		scalar coNumRatio = normPrev/normInit;

		coNumRatio = max(min(coNumRatio, pseudoCoNumMaxIncr), pseudoCoNumMinDecr);

		if (localTimestepping)
		{
			pseudoCoField *= coNumRatio;
			pseudoCoField = max(min(pseudoCoField, pseudoCoNumMax), pseudoCoNumMin);
			pseudoCoField.correctBoundaryConditions();
		}
		else
		{
			pseudoCoNum *= coNumRatio;
			pseudoCoNum.value() = max(min(pseudoCoNum.value(), pseudoCoNumMax), pseudoCoNumMin);
		}
	}

	prevRes.clear();
	prevRes.set(new residualsIO(initRes()));
}

tmp< volScalarField > gamma = pThermo->gamma();
const volScalarField& psi = pThermo->psi();
tmp< surfaceScalarField > lambda;

if (mesh.moving())
{
	lambda = fvc::interpolate(sqrt(gamma()/psi))
				+ mag((fvc::interpolate(U) & mesh.Sf()/mesh.magSf())
				- flux.MRFFaceVelocity()
				- fvc::meshPhi(U)/mesh.magSf());
}
else
{
	 lambda = fvc::interpolate(sqrt(gamma()/psi))
				+ mag((fvc::interpolate(U) & mesh.Sf()/mesh.magSf())
				- flux.MRFFaceVelocity());
}

if (localTimestepping)
{
	rPseudoDeltaT == dimensionedScalar("0", dimless/dimTime, 0.0);

	surfaceScalarField frdt = mesh.nonOrthDeltaCoeffs()*(lambda());

	const labelList& own = mesh.owner();
	const labelList& nei = mesh.neighbour();

	forAll(own,facei)
	{
		rPseudoDeltaT[own[facei]] = max(rPseudoDeltaT[own[facei]],frdt[facei]);
		rPseudoDeltaT[nei[facei]] = max(rPseudoDeltaT[nei[facei]],frdt[facei]);
	}

	forAll(frdt.boundaryField(), patchi)
	{
		const labelUList& fc = mesh.boundary()[patchi].faceCells();

		if (mesh.boundary()[patchi].coupled())
		{
			forAll(fc, bfacei)
			{
				rPseudoDeltaT[fc[bfacei]] =
					max
					(
						rPseudoDeltaT[fc[bfacei]],
						frdt.boundaryField()[patchi][bfacei]
					);
			}
		}
		else if
		(
			mesh.boundary()[patchi].type() == "wall"
		)
		{
			const scalarField& pNonOrthDeltaCoef(mesh.nonOrthDeltaCoeffs().boundaryField()[patchi]);

			const vectorField& pSf(mesh.Sf().boundaryField()[patchi]);
			const scalarField& pMagSf(mesh.magSf().boundaryField()[patchi]);
			const scalarField& pMRFFaceVel(flux.MRFFaceVelocity().boundaryField()[patchi]);

			scalarField pMeshPhi(pMagSf.size());

			if (mesh.moving())
			{
				pMeshPhi = fvc::meshPhi(U)().boundaryField()[patchi];
			}

			forAll(fc, bfacei)
			{
				scalar pLambda;

				if (mesh.moving())
				{
					pLambda =
							0.5*pNonOrthDeltaCoef[bfacei]*
							(Foam::sqrt(gamma()[fc[bfacei]]/psi[fc[bfacei]])
								+ mag((U[fc[bfacei]] & pSf[bfacei]/pMagSf[bfacei])
								- pMRFFaceVel[bfacei]
								- pMeshPhi[bfacei]/pMagSf[bfacei]));
				}
				else
				{
					 pLambda =
							0.5*pNonOrthDeltaCoef[bfacei]*
							 (Foam::sqrt(gamma()[fc[bfacei]]/psi[fc[bfacei]])
					 	 	 	 + mag((U[fc[bfacei]] & pSf[bfacei]/pMagSf[bfacei])
								 - pMRFFaceVel[bfacei]));
				}



				rPseudoDeltaT[fc[bfacei]] =
					max
					(
						rPseudoDeltaT[fc[bfacei]], // Cell value based on max of internal neighbours
						pLambda
					);
			}
		}
	}

	rPseudoDeltaT /= pseudoCoField;
	rPseudoDeltaT.correctBoundaryConditions();
}
else
{
	rPseudoDeltaT = max(mesh.deltaCoeffs()*lambda); // Global timestep
	rPseudoDeltaT /= pseudoCoNum;
}

if (localTimestepping)
{
	scalar totCells = mesh.globalData().nTotalCells();

	Info<< "Pseudo Courant No: "
		<< "Min: "
		<< gMin(pseudoCoField) 
		<< " Mean: "
		<< gSum(pseudoCoField)/totCells
		<< " Max: "
		<< gMax(pseudoCoField) << endl;

	Info<< "Pseudo deltaT: "
		<< "Min: "
		<< 1.0/gMax(rPseudoDeltaT)
		<< " Mean: "
		<< gSum((1.0/rPseudoDeltaT)())/totCells
		<< " Max: "
		<< 1.0/gMin(rPseudoDeltaT) << endl;
}
else
{
	Info << "Pseudo Courant No: " << pseudoCoNum.value() << endl;
	Info << "Pseudo deltaT: " << 1.0/gMin(rPseudoDeltaT) << endl;
}
