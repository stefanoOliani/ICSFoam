   // Detect steady-state analysis
    const dictionary& ddtControls =
    		mesh.schemesDict().subDict("ddtSchemes");

    wordList ddtToc (ddtControls.toc());

    bool steadyState = false;

    forAll(ddtToc,s)
    {
        word ddtScheme(ddtToc[s]);
        word ddtSchemeLastWord;
        const tokenList& tokens = ddtControls.lookup(ddtToc[s]);
        if (tokens.last().isWord() && tokens.last().wordToken() == "steadyState")
        {
            if (ddtToc[s] == "default" || ddtToc[s] == "rhoU")
            {
                steadyState = true;
            }
        }
    }

    if (steadyState)
    {
        Info << "Steady-state analysis detected" << nl << endl;
    }
    else
    {
        Info << "Transient analysis detected" << nl << endl;
    }

    pseudotimeControl solnControl(mesh, steadyState, 2, 1);

    if (steadyState)
    {
        solnControl.setCorr(runTime.startTimeIndex());
    }

    bool localTimestepping = solnControl.dict().lookupOrDefault<Switch>
							(
								"localTimestepping",
								true
							);

    bool localTimesteppingBounding(false);
    scalar localTimesteppingLowerBound = 0.0;

    if (localTimestepping)
    {
        Info << "Local timestepping selected" << nl << endl;

        localTimesteppingBounding =
            solnControl.dict().lookupOrDefault<Switch>
            (
                "localTimesteppingBounding",
                true
            );

        localTimesteppingLowerBound =
            solnControl.dict().lookupOrDefault<scalar>
            (
                "localTimesteppingLowerBound",
                0.95
            );

        localTimesteppingLowerBound =
            (localTimesteppingLowerBound > 0 ? localTimesteppingLowerBound : 0.0);

        localTimesteppingLowerBound =
            (localTimesteppingLowerBound < 0.99 ? localTimesteppingLowerBound : 0.99);
    }
    else
    {
        Info << "Global timestepping selected" << nl << endl;
    }

    #include "createFields.H"

    // Read or initialise pseudo Co number
    if (!localTimestepping)
    {
        Info << "Initial pseudo Courant No: " << pseudoCoNum.value() << nl << endl;
    }
    else
    {
        scalar totCells = mesh.globalData().nTotalCells();

        Info<< "Initial pseudo Courant No: "
            << "Min: "
            << min(pseudoCoField).value()
            << " Mean: "
            << sum(pseudoCoField).value()/totCells
            << " Max: "
            << max(pseudoCoField).value() << nl << endl;
    }
