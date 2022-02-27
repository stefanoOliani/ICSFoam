// Set all the meshes to the correct position to
// compute the grid fluxes

scalar timevalue_temp = runTime.value();
label index_temp = runTime.timeIndex();
dimensionedScalar deltaT_temp = runTime.deltaT();

forAll(subTimeMeshes,K)
{
	  runTime.setTime(K*deltaPeriod + deltaPeriod/5000, index_temp + 1);
	  runTime.setDeltaT(deltaPeriod/5000, false);

	  subTimeMeshes[K].update();
}

runTime.setTime(timevalue_temp,index_temp);
runTime.setDeltaT(deltaT_temp);
