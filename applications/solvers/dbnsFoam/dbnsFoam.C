/*---------------------------------------------------------------------------*\

    ICSFoam: a library for Implicit Coupled Simulations in OpenFOAM
  
    Copyright (C) 2022  Stefano Oliani

    https://turbofe.it

-------------------------------------------------------------------------------
License
    This file is part of ICSFOAM.

    ICSFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFOAM.  If not, see <http://www.gnu.org/licenses/>.


Author
    Stefano Oliani
    Fluid Machinery Research Group, University of Ferrara, Italy
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "pseudotimeControl.H"
#include "dynamicFvMesh.H"
#include "IOMRFCoupledZoneList.H"
#include "IOMRFTranslatingZoneList.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "fvcSmooth.H"
#include "coupledMatrix.H"
#include "convectiveFluxScheme.H"
#include "viscousFluxScheme.H"
#include "fvOptions.H"

#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient db solver for compressible turbulent flow.\n"
    );


	#include "addCheckCaseOptions.H"
	#include "setRootCaseLists.H"
	#include "createTime.H"
	#include "createDynamicFvMesh.H"

	#include "initialise.H"

    if (!inviscid)
    {
        turbulence->validate();
    }

    if (!steadyState)
    {
        #include "createTimeControls.H"

        #include "compressibleCourantNo.H"

        CoNum = maxCo/(maxCo/stabilise(CoNum, SMALL));

        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (steadyState)
    {
        Info << "\nStarting pseudotime iteration loop\n" << endl;
    }
    else
    {
        Info<< "\nStarting time loop\n" << endl;
    }

    while (runTime.run())
    {
        if (!steadyState)
        {
            #include "createTimeControls.H"

            #include "compressibleCourantNo.H"

            CoNum = maxCo/(maxCo/stabilise(CoNum, SMALL));

            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

	#include "beginTimeStep.H"

        // --- Outer corrector loop
        bool notFinished;

        while ( (notFinished = solnControl.loop()) )
        {
	    #include "outerLoop.H"

            // Correct turbulence
            turbulence->correct();

            // Update fields
            #include "updateFields.H"

            if (steadyState)
            {
                break;
            }
        }

        if (scalarTransport)
	{
            int nIterScalarTransp = scalarDict.get<int>("nIterations");

            for (int i = 0;  i<nIterScalarTransp ; i++)
            {
		#include "calculateScalarTransport.H"
            }
	}

        volScalarField rhoRes(fvc::div(phi));
	scalar L2NormRho = Foam::sqrt(sum(sqr(rhoRes.primitiveField()*mesh.V())));
	Info<< "rho L2 Residual: "<< Foam::log10(L2NormRho)  << endl;

        if (steadyState && !notFinished)
        {
            runTime.writeAndEnd();
        }
        else
        {
            runTime.write();
        }

        runTime.printExecutionTime(Info);

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
