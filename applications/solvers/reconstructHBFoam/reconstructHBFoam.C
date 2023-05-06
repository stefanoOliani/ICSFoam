/*---------------------------------------------------------------------------*\

    Copyright (C) 2004-2011 OpenCFD Ltd.
    Copyright (C) 2022 Stefano Oliani

-------------------------------------------------------------------------------
License
    This file is part of ICSFoam.

    ICSFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    ICSFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with ICSFoam.  If not, see <http://www.gnu.org/licenses/>.

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

#include "complex.H"
#include "IOHBZoneList.H"

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
	#include "createMeshes.H"

	#include "createFields.H"

    if (!inviscid)
    {
    	forAll(subTimeMeshes,K)
    	{
			turbulenceHB[K].validate();
    	}
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting iteration loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mesh.update();

        scalar actTime = runTime.value();

        HB.reconstruct(U, UHB, actTime, true);
        HB.reconstruct(p, pHB, actTime);
        HB.reconstruct(T, THB, actTime);
        HB.reconstruct(rho, rhoHB, actTime);

        p.correctBoundaryConditions();
        U.correctBoundaryConditions();
        T.correctBoundaryConditions();

        phi = mesh.Sf() & fvc::interpolate(rho*U);

        if (mesh.moving())
        {
            URel = U - fvc::reconstruct(fvc::meshPhi(U));
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
