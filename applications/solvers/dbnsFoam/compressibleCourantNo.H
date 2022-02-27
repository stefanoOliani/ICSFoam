scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (mesh.nInternalFaces())
{
    tmp< surfaceScalarField > lambda;

    if (mesh.moving())
    {
        lambda = fvc::interpolate(sqrt(thermo.gamma()/thermo.psi()))
        		+ mag((fvc::interpolate(U)&mesh.Sf()/mesh.magSf())
        		- flux.MRFFaceVelocity()
        		- fvc::meshPhi(U)/mesh.magSf());
    }
    else
    {
        lambda = fvc::interpolate(sqrt(thermo.gamma()/thermo.psi()))
        		+ mag((fvc::interpolate(U)&mesh.Sf()/mesh.magSf())
        		- flux.MRFFaceVelocity());
    }

    surfaceScalarField amaxSfbyDelta
    (
        lambda*mesh.deltaCoeffs() 
    );

    CoNum = gMax(amaxSfbyDelta)*runTime.deltaTValue();

    meanCoNum =
            (gSum((amaxSfbyDelta*mesh.magSf())())/gSum(mesh.magSf()))
           *runTime.deltaTValue();
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << " " << CoNum << endl;

// ************************************************************************* //
