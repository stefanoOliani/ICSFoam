    const dictionary& dict = subTimeMeshes[0].solutionDict().subDict("pseudoTime");

	volScalarField& p_ = pThermo[K].p();
	volScalarField& e_ = pThermo[K].he();
	volScalarField& T_ = pThermo[K].T();
	const volScalarField& psi_ = pThermo[K].psi();

	volScalarField& rhoIncr = scalarSolIncrement[2*K];
	volVectorField& rhoUIncr = vectorSolIncrement[K];
	volScalarField& rhoEIncr = scalarSolIncrement[2*K+1];

	rho_.ref() += rhoIncr.internalField();
	rhoU_.ref() += rhoUIncr.internalField();
	rhoE_.ref() += rhoEIncr.internalField();

    dimensionedScalar rhoMin
    (
        "rhoMin", 
        dimDensity, 
        dict.lookupOrDefault<scalar>("rhoMin", -GREAT)
    );
    bound(rho_, rhoMin);

    // Compute internal fields

    U_.ref() = rhoU_.internalField()/rho_.internalField();

    e_.ref() = rhoE_.internalField()/rho_.internalField()
              - 0.5*magSqr(U_.internalField());

    dimensionedScalar TMin
    (
        "TMin", 
        dimTemperature, 
        dict.lookupOrDefault<scalar>("TMin", SMALL)
    );
    dimensionedScalar TMax
    (
        "TMax", 
        dimTemperature, 
        dict.lookupOrDefault<scalar>("TMax", GREAT)
    );
    
    volScalarField TBound
    (
        IOobject("TBound", runTime.timeName(), subTimeMeshes[K]),
        subTimeMeshes[K],
        TMin
    );
    volScalarField eBound = pThermo[K].he(p_, TBound);

    if (max(neg(e_-eBound)).value() > 0.5)
    {
        e_ = max(e_, eBound);

         Info<< "Bounding " << e_.name()
             << " to TMin: " << TMin.value()
             << endl;
    }

    if (TMax.value() < GREAT)
    {
        TBound == TMax;
        eBound == pThermo[K].he(p_, TBound);
        if (max(pos(e_-eBound)).value() > 0.5)
        {
            e_ = min(e_, eBound);
            Info<< "Bounding " << e_.name()
                << " to TMax: " << TMax.value()
                << endl;
        }
    }

    pThermo[K].correct();

    p_.ref() = rho_.internalField()/psi_.internalField();

    rhoU_.ref() = rho_.internalField()*U_.internalField();

    rhoE_.ref() =
        rho_.internalField()*
        (
            e_.internalField() + 0.5*magSqr(U_.internalField())
        );

    p_.correctBoundaryConditions();         
    U_.correctBoundaryConditions();
    MRF[K].correctBoundaryVelocity(U_);
    T_.correctBoundaryConditions();

    T_.boundaryFieldRef() == max(T_.boundaryField(), TMin.value());

    if (TMax.value() < GREAT)
    {
        T_.boundaryFieldRef() == min(T_.boundaryField(), TMax.value());
    }

    e_.boundaryFieldRef() == pThermo[K].he(p_, T_)->boundaryField();

    pThermo[K].correct();                     

    rho_.boundaryFieldRef() == psi_.boundaryField()*p_.boundaryField();

    rhoU_.boundaryFieldRef() = rho_.boundaryField()*U_.boundaryField();

    rhoE_.boundaryFieldRef() =
        rho_.boundaryField()*
        (
            e_.boundaryField()
          + 0.5*magSqr(U_.boundaryField())
        );

