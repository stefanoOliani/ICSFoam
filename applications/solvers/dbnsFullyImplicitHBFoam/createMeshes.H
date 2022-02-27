// Read the number of instants from fvSolution
fvSolution solutionDict(runTime);

// Read harmonics coeffs dict
dictionary harmonics = solutionDict.subDict("harmonicBalance");

int nO = harmonics.getOrDefault<label>("instantsNumber",3);

PtrList<dynamicFvMesh> subTimeMeshes(nO);

if (!Pstream::parRun())
{
	if(fileHandler().exists(cwd()/"tmpCnst"))
	{
		rmDir(cwd()/"tmpCnst");
	}

	if(fileHandler().exists(cwd()/"tmpSys"))
	{
		rmDir(cwd()/"tmpSys");
	}

	forAll(subTimeMeshes, i)
	{
		fileName subLevConst = runTime.constant()/"subTimeLevel" + Foam::name(i);

		if(fileHandler().exists(subLevConst))
		{
			rmDir(subLevConst);
		}

		fileName subLevSys = runTime.system()/"subTimeLevel" + Foam::name(i);

		if(fileHandler().exists(subLevSys))
		{
			rmDir(subLevSys);
		}
	}

	cp(runTime.constant(), "tmpCnst");
	cp(runTime.system(), "tmpSys");

	forAll(subTimeMeshes, i)
	{
		word itrName = Foam::name(i);

		Info<< " \n Create mesh for sub-time level " << itrName << nl << endl;

		cp("tmpCnst",runTime.constant()/"subTimeLevel" + itrName);

		cp("tmpSys",runTime.system()/"subTimeLevel" + itrName);

		subTimeMeshes.set
		(
		  i,
			autoPtr<dynamicFvMesh>
			(
				dynamicFvMesh::New
				(
					args,
					runTime
				)
			)
		);

		subTimeMeshes[i].polyMesh::rename("subTimeLevel" + itrName);
	}
}
else
{
	forAll(subTimeMeshes, i)
	{
		word itrName = Foam::name(i);

		Info<< " \n Create mesh for sub-time level " << itrName << nl << endl;

		subTimeMeshes.set
		(
		  i,
			autoPtr<dynamicFvMesh>
			(
				dynamicFvMesh::New
				(
					args,
					runTime
				)
			)
		);

		subTimeMeshes[i].polyMesh::rename("subTimeLevel" + itrName);
	}
}

IOHBZoneList HB(subTimeMeshes[0]);
HB.setInstants();

scalar Time_temp = runTime.value();

const scalarList& snapshots = HB.selectedSnapshots();

const scalar deltaPeriod = snapshots[1] - snapshots[0];

forAll(subTimeMeshes, i)
{
	runTime.value()=i*deltaPeriod;
	subTimeMeshes[i].update();
}

runTime.value() = Time_temp;

if (!Pstream::parRun())
{
	rmDir(cwd()/"tmpCnst");
	rmDir(cwd()/"tmpSys");
}

