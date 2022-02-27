    // Clear out residuals from previous time step
    if (!steadyState)
    {
    	initRes.clear();
        prevRes.clear();
    }

    	bool moveMeshOuterCorrectors =
            solnControl.dict().lookupOrDefault<Switch>
            (
                "moveMeshOuterCorrectors",
                false
            );

        // Reset pseudo Co number before every time step
        if (solnControl.dict().lookupOrDefault<bool>("resetPseudo", false) )
        {
            if (!localTimestepping)
            {
                pseudoCoNum.value() =
                    solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
            }
            else
            {
                pseudoCoField == solnControl.dict().lookupOrDefault<scalar>("pseudoCoNum", 1);
            }
        }

        scalar pseudoCoNumMin =
            solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMin", 0.1);

        scalar pseudoCoNumMax =
            solnControl.dict().lookupOrDefault<scalar>("pseudoCoNumMax", 25);

        scalar pseudoCoNumMaxIncr = solnControl.dict().lookupOrDefault<scalar>
        (
            "pseudoCoNumMaxIncreaseFactor",
            1.25
        );

        scalar pseudoCoNumMinDecr = solnControl.dict().lookupOrDefault<scalar>
        (
            "pseudoCoNumMinDecreaseFactor",
            0.1
        );
