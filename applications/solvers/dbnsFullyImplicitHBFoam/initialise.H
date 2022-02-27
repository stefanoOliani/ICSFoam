    Info << "Harmonic Balance analysis detected" << nl << endl;

    PtrList<pseudotimeControl> solnControl(nO);

    forAll(solnControl,i)
    {
    	solnControl.set
    	(
    		i,
    	    new pseudotimeControl
    	    (
    			subTimeMeshes[i],
				true,
				2,
				1
    	    )
    	);

        solnControl[i].setCorr(runTime.startTimeIndex());
    }

    bool localTimesteppingBounding(false);
    scalar localTimesteppingLowerBound = 0.0;

	Info << "Local timestepping selected" << nl << endl;

	localTimesteppingBounding =
		solnControl[0].dict().lookupOrDefault<Switch>
		(
			"localTimesteppingBounding",
			true
		);

	localTimesteppingLowerBound =
		solnControl[0].dict().lookupOrDefault<scalar>
		(
			"localTimesteppingLowerBound",
			0.95
		);

	localTimesteppingLowerBound =
		(localTimesteppingLowerBound > 0 ? localTimesteppingLowerBound : 0.0);

	localTimesteppingLowerBound =
		(localTimesteppingLowerBound < 0.99 ? localTimesteppingLowerBound : 0.99);

    #include "createFields.H"

	forAll(subTimeMeshes,i)
	{
	    // Read or initialise pseudo Co number
	    // NOTE: It is not necessarily read before every outer iteration
	    // (see resetPseudo)
		scalar totCells = subTimeMeshes[i].globalData().nTotalCells();

		Info<<"Sub-time mesh "<<i<< nl
			<< "Initial pseudo Courant No: "
			<< "Min: "
			<< min(pseudoCoField[i]).value()
			<< " Mean: "
			<< sum(pseudoCoField[i]).value()/totCells
			<< " Max: "
			<< max(pseudoCoField[i]).value() << nl << endl;
	}

