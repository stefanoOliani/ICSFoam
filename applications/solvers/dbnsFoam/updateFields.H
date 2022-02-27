// Optional failsafe to avoid -ve rho. Avoid dividing by very small quantity
const dictionary& dict = mesh.solutionDict().subDict("pseudoTime");

volScalarField& rhoIncr = scalarSolIncrement[0];
    volVectorField& rhoUIncr = vectorSolIncrement[0];
    volScalarField& rhoEIncr = scalarSolIncrement[1];

    rho.ref() += rhoIncr.internalField();
    rhoU.ref() += rhoUIncr.internalField();
    rhoE.ref() += rhoEIncr.internalField();

dimensionedScalar rhoMin
(
    "rhoMin", 
    dimDensity, 
    dict.lookupOrDefault<scalar>("rhoMin", -GREAT)
);
bound(rho, rhoMin);

// Compute internal fields

U.ref() = rhoU.internalField()/rho.internalField();

e.ref() = rhoE.internalField()/rho.internalField()
          - 0.5*magSqr(U.internalField());

// Temperature bound
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

// Bound energy
volScalarField TBound
(
    IOobject("TBound", runTime.timeName(), mesh), 
    mesh,
    TMin
);
volScalarField eBound = thermo.he(thermo.p(), TBound);

if (max(neg(e-eBound)).value() > 0.5)
{
    e = max(e, eBound);

     Info<< "Bounding " << e.name()
         << " to TMin: " << TMin.value()
         << endl;
}

// Only do max bound if it was specified
if (TMax.value() < GREAT)
{
    TBound == TMax;
    eBound == thermo.he(thermo.p(), TBound);
    if (max(pos(e-eBound)).value() > 0.5)
    {
        e = min(e, eBound);
        Info<< "Bounding " << e.name()
            << " to TMax: " << TMax.value()
            << endl;
    }
}

// Calc T and psi from e
thermo.correct();

p.ref() = rho.internalField()/psi.internalField();

// Recalc rhoU and rhoE in case rho or e were bounded
rhoU.ref() = rho.internalField()*U.internalField();

rhoE.ref() =
    rho.internalField()*
    (
		e.internalField() + 0.5*magSqr(U.internalField())
    );

// Correct boundary fields
p.correctBoundaryConditions();          
U.correctBoundaryConditions();
T.correctBoundaryConditions();

T.boundaryFieldRef() == max(T.boundaryField(), TMin.value());

if (TMax.value() < GREAT)
{
    T.boundaryFieldRef() == max(T.boundaryField(), TMax.value());
}

e.boundaryFieldRef() == thermo.he(p, T)->boundaryField();

thermo.correct();                      

rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField(); 
rho.correctBoundaryConditions();        

rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();

rhoE.boundaryFieldRef() =
    rho.boundaryField()*
    (
        e.boundaryField()
      + 0.5*magSqr(U.boundaryField())
    );

