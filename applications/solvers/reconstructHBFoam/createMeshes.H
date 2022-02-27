// Read the number of instants from fvSolution
fvSolution solutionDict(runTime);

// Read harmonics coeffs dict
dictionary harmonics = solutionDict.subDict("harmonicBalance");

int nO = harmonics.getOrDefault<label>("instantsNumber",3);


rmDir("tmpCnst");
rmDir("tmpSys");

PtrList<dynamicFvMesh> subTimeMeshes(nO);

forAll(subTimeMeshes, i)
{
	word itrName = Foam::name(i);
	rmDir(runTime.constant()/"subTimeLevel" + itrName);
	rmDir(runTime.system()/"subTimeLevel" + itrName);
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

rmDir("tmpCnst");
rmDir("tmpSys");

