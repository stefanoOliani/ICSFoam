volVectorField& U_ = U[K];

// Switched Evolution Relaxation
if (initRes.valid())
{
	if (!solnControl[K].firstIter() && prevRes.valid())
	{
		pseudoCoField[K] *= coNumRatio;

		pseudoCoField[K] = max(min(pseudoCoField[K], pseudoCoNumMax), pseudoCoNumMin);
		pseudoCoField[K].correctBoundaryConditions();
	}
}

tmp< volScalarField > gamma = pThermo[K].gamma();
const volScalarField& psi = pThermo[K].psi();
tmp< surfaceScalarField > lambda;

if (mesh.moving())
{
	lambda = fvc::interpolate(sqrt(gamma()/psi))
				+ mag((fvc::interpolate(U_) & mesh.Sf()/mesh.magSf())
				- convectiveFlux[K].MRFFaceVelocity()
				- fvc::meshPhi(U_)/mesh.magSf());
}
else
{
	 lambda = fvc::interpolate(sqrt(gamma()/psi))
				+ mag((fvc::interpolate(U_) & mesh.Sf()/mesh.magSf())
				- convectiveFlux[K].MRFFaceVelocity());
}

rPseudoDeltaT[K] == dimensionedScalar("0", dimless/dimTime, 0.0);

surfaceScalarField frdt = mesh.nonOrthDeltaCoeffs()*(lambda());

const labelList& own = mesh.owner();
const labelList& nei = mesh.neighbour();

forAll(own,facei)
{
	rPseudoDeltaT[K][own[facei]] = max(rPseudoDeltaT[K][own[facei]],frdt[facei]);
	rPseudoDeltaT[K][nei[facei]] = max(rPseudoDeltaT[K][nei[facei]],frdt[facei]);
}

forAll(frdt.boundaryField(), patchi)
{
	const labelUList& fc = mesh.boundary()[patchi].faceCells();

	if (mesh.boundary()[patchi].coupled())
	{
		forAll(fc, bfacei)
		{
			rPseudoDeltaT[K][fc[bfacei]] =
				max
				(
					rPseudoDeltaT[K][fc[bfacei]],
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
		const scalarField& pMRFFaceVel(convectiveFlux[K].MRFFaceVelocity().boundaryField()[patchi]);

		forAll(fc, bfacei)
		{
			scalar pLambda;

			if (mesh.moving())
			{
				const scalarField& pMeshPhi(fvc::meshPhi(U_)().boundaryField()[patchi]);

				pLambda =
						0.5*pNonOrthDeltaCoef[bfacei]*
						(Foam::sqrt(gamma()[fc[bfacei]]/psi[fc[bfacei]])
							+ mag((U_[fc[bfacei]] & pSf[bfacei]/pMagSf[bfacei])
							- pMRFFaceVel[bfacei]
							- pMeshPhi[bfacei]/pMagSf[bfacei]));
			}
			else
			{
				 pLambda =
						0.5*pNonOrthDeltaCoef[bfacei]*
						 (Foam::sqrt(gamma()[fc[bfacei]]/psi[fc[bfacei]])
							 + mag((U_[fc[bfacei]] & pSf[bfacei]/pMagSf[bfacei])
							 - pMRFFaceVel[bfacei]));
			}



			rPseudoDeltaT[K][fc[bfacei]] =
				max
				(
					rPseudoDeltaT[K][fc[bfacei]], // Cell value based on max of internal neighbours
					pLambda
				);
		}
	}
}

rPseudoDeltaT[K] /= pseudoCoField[K];
rPseudoDeltaT[K].correctBoundaryConditions();

scalar totCells = mesh.globalData().nTotalCells();

Info<< " \n Time level "<< K << "\n"<<endl;

Info<< "Pseudo Courant No: "
	<< "Min: "
	<< min(pseudoCoField[K]).value()
	<< " Mean: "
	<< sum(pseudoCoField[K]).value()/totCells
	<< " Max: "
	<< max(pseudoCoField[K]).value() << endl;

Info<< "Pseudo deltaT: "
	<< "Min: "
	<< 1.0/max(rPseudoDeltaT[K]).value()
	<< " Mean: "
	<< sum(1.0/rPseudoDeltaT[K]).value()/totCells
	<< " Max: "
	<< 1.0/min(rPseudoDeltaT[K]).value() << endl;

