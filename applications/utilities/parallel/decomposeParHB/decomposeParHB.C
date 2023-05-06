/*---------------------------------------------------------------------------*\

    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd
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

Application
    decomposeParHB

Group
    grpParallelUtilities

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposeParHB [OPTIONS]

    Options:
      - \par -allRegions
        Decompose all regions in regionProperties. Does not check for
        existence of processor*.

      - \par -case \<dir\>
        Specify case directory to use (instead of the cwd).

      - \par -cellDist
        Write the cell distribution as a labelList, for use with 'manual'
        decomposition method and as a volScalarField for visualization.

      - \par -constant
        Include the 'constant/' dir in the times list.

      - \par -copyUniform
        Copy any \a uniform directories too.

      - \par -copyZero
        Copy \a 0 directory to processor* rather than decompose the fields.

      - \par -debug-switch \<name=val\>
        Specify the value of a registered debug switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -decomposeParHBDict \<file\>
        Use specified file for decomposeParHB dictionary.

      - \par -dry-run
        Test without writing the decomposition. Changes -cellDist to
        only write volScalarField.

      - \par -fields
        Use existing geometry decomposition and convert fields only.

      - \par fileHandler \<handler\>
        Override the file handler type.

      - \par -force
        Remove any existing \a processor subdirectories before decomposing the
        geometry.

      - \par -ifRequired
        Only decompose the geometry if the number of domains has changed from a
        previous decomposition. No \a processor subdirectories will be removed
        unless the \a -force option is also specified. This option can be used
        to avoid redundant geometry decomposition (eg, in scripts), but should
        be used with caution when the underlying (serial) geometry or the
        decomposition method etc. have been changed between decompositions.

      - \par -info-switch \<name=val\>
        Specify the value of a registered info switch. Default is 1
        if the value is omitted. (Can be used multiple times)

      - \par -latestTime
        Select the latest time.

      - \par -lib \<name\>
        Additional library or library list to load (can be used multiple times).

      - \par -noFunctionObjects
        Do not execute function objects.

      - \par -noSets
        Skip decomposing cellSets, faceSets, pointSets.

      - \par -noZero
        Exclude the \a 0 dir from the times list.

      - \par -opt-switch \<name=val\>
        Specify the value of a registered optimisation switch (int/bool).
        Default is 1 if the value is omitted. (Can be used multiple times)

      - \par -region \<regionName\>
        Decompose named region. Does not check for existence of processor*.

      - \par -time \<ranges\>
        Override controlDict settings and decompose selected times. Does not
        re-decompose the mesh i.e. does not handle moving mesh or changing
        mesh cases. Eg, ':10,20 40:70 1000:', 'none', etc.

      - \par -verbose
        Additional verbosity.

      - \par -doc
        Display documentation in browser.

      - \par -doc-source
        Display source code in browser.

      - \par -help
        Display short help and exit.

      - \par -help-man
        Display full help (manpage format) and exit.

      - \par -help-notes
        Display help notes (description) and exit.

      - \par -help-full
        Display full help and exit.

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "IOobjectList.H"
#include "domainDecomposition.H"
#include "labelIOField.H"
#include "labelFieldIOField.H"
#include "scalarIOField.H"
#include "scalarFieldIOField.H"
#include "vectorIOField.H"
#include "vectorFieldIOField.H"
#include "sphericalTensorIOField.H"
#include "sphericalTensorFieldIOField.H"
#include "symmTensorIOField.H"
#include "symmTensorFieldIOField.H"
#include "tensorIOField.H"
#include "tensorFieldIOField.H"
#include "pointFields.H"
#include "regionProperties.H"

#include "readFields.H"
#include "dimFieldDecomposer.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "decompositionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

const labelIOList& procAddressing
(
    const PtrList<fvMesh>& procMeshList,
    const label proci,
    const word& name,
    PtrList<labelIOList>& procAddressingList
)
{
    const fvMesh& procMesh = procMeshList[proci];

    if (!procAddressingList.set(proci))
    {
        procAddressingList.set
        (
            proci,
            new labelIOList
            (
                IOobject
                (
                    name,
                    procMesh.facesInstance(),
                    procMesh.meshSubDir,
                    procMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
    }
    return procAddressingList[proci];
}


void decomposeUniform
(
    const bool copyUniform,
    const domainDecomposition& mesh,
    const Time& processorDb,
    const word& regionDir = word::null
)
{
    const Time& runTime = mesh.time();

    // Any uniform data to copy/link?
    const fileName uniformDir(regionDir/"uniform");

    if (fileHandler().isDir(runTime.timePath()/uniformDir))
    {
        Info<< "Detected additional non-decomposed files in "
            << runTime.timePath()/uniformDir
            << endl;

        const fileName timePath =
            fileHandler().filePath(processorDb.timePath());

        // If no fields have been decomposed the destination
        // directory will not have been created so make sure.
        mkDir(timePath);

        if (copyUniform || mesh.distributed())
        {
            if (!fileHandler().exists(timePath/uniformDir))
            {
                fileHandler().cp
                (
                    runTime.timePath()/uniformDir,
                    timePath/uniformDir
                );
            }
        }
        else
        {
            // Link with relative paths
            string parentPath = string("..")/"..";

            if (regionDir != word::null)
            {
                parentPath = parentPath/"..";
            }

            fileName currentDir(cwd());
            chDir(timePath);

            if (!fileHandler().exists(uniformDir))
            {
                fileHandler().ln
                (
                    parentPath/runTime.timeName()/uniformDir,
                    uniformDir
                );
            }
            chDir(currentDir);
        }
    }
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Decompose a mesh and fields of a case for parallel execution"
    );

    argList::noParallel();
    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Use specified file for decomposePar dictionary"
    );
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "Operate on all regions in regionProperties"
    );
    argList::addBoolOption
    (
        "dry-run",
        "Test without writing the decomposition. "
        "Changes -cellDist to only write volScalarField."
    );
    argList::addBoolOption
    (
        "verbose",
        "Additional verbosity"
    );
    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method and as a volScalarField for visualization."
    );
    argList::addBoolOption
    (
        "copyZero",
        "Copy 0/ directory to processor*/ rather than decompose the fields"
    );
    argList::addBoolOption
    (
        "copyUniform",
        "Copy any uniform/ directories too"
    );
    argList::addBoolOption
    (
        "fields",
        "Use existing geometry decomposition and convert fields only"
    );
    argList::addBoolOption
    (
        "noSets",
        "Skip decomposing cellSets, faceSets, pointSets"
    );
    argList::addBoolOption
    (
        "force",
        "Remove existing processor*/ subdirs before decomposing the geometry"
    );
    argList::addBoolOption
    (
        "ifRequired",
        "Only decompose geometry if the number of domains has changed"
    );

    // Allow explicit -constant, have zero from time range
    timeSelector::addOptions(true, false);  // constant(true), zero(false)

    #include "setRootCase.H"

    const bool dryrun           = args.found("dry-run");
//    const bool optRegion        = args.found("region");
//    const bool allRegions       = args.found("allRegions");
    const bool writeCellDist    = args.found("cellDist");
 //   const bool verbose          = args.found("verbose");

    // Most of these are ignored for dry-run (not triggered anywhere)
    const bool copyZero         = args.found("copyZero");
    const bool copyUniform      = args.found("copyUniform");
    const bool decomposeSets    = !args.found("noSets");
    const bool decomposeIfRequired = args.found("ifRequired");

    bool decomposeFieldsOnly = args.found("fields");
    bool forceOverwrite      = args.found("force");


    // Set time from database
    #include "createTime.H"

    // Allow override of time (unless dry-run)
    instantList times;
    if (dryrun)
    {
        Info<< "\ndry-run: ignoring -copy*, -fields, -force, time selection"
            << nl;
    }
    else
    {
        times = timeSelector::selectIfPresent(runTime, args);
    }

    // Allow override of decomposeParHBDict location
    fileName decompDictFile(args.get<fileName>("decomposeParDict", ""));
    if (!decompDictFile.empty() && !decompDictFile.isAbsolute())
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    Info<< "\n\nDecomposing base mesh " << nl << endl;

	// Determine the existing processor count directly
	label nProcs = fileHandler().nProcs(runTime.path(), word::null);

	// Get requested numberOfSubdomains directly from the dictionary.
	// Note: have no mesh yet so cannot use decompositionModel::New
	const label nDomains = decompositionMethod::nDomains
	(
		IOdictionary
		(
			IOobject::selectIO
			(
				IOobject
				(
					decompositionModel::canonicalName,
					runTime.time().system(),
					word::null,  // region (if non-default)
					runTime,
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				),
				decompDictFile
			)
		)
	);

	// Give file handler a chance to determine the output directory
	const_cast<fileOperation&>(fileHandler()).setNProcs(nDomains);

	if (decomposeFieldsOnly)
	{
		// Sanity check on previously decomposed case
		if (nProcs != nDomains)
		{
			FatalErrorInFunction
				<< "Specified -fields, but the case was decomposed with "
				<< nProcs << " domains"
				<< nl
				<< "instead of " << nDomains
				<< " domains as specified in decomposeParHBDict" << nl
				<< exit(FatalError);
		}
	}
	else if (nProcs)
	{
		bool procDirsProblem = true;

		if (decomposeIfRequired && nProcs == nDomains)
		{
			// We can reuse the decomposition
			decomposeFieldsOnly = true;
			procDirsProblem = false;
			forceOverwrite = false;

			Info<< "Using existing processor directories" << nl;
		}

		if (forceOverwrite)
		{
			Info<< "Removing " << nProcs
				<< " existing processor directories" << endl;

			// Remove existing processors directory
			fileNameList dirs
			(
				fileHandler().readDir
				(
					runTime.path(),
					fileName::Type::DIRECTORY
				)
			);
			forAllReverse(dirs, diri)
			{
				const fileName& d = dirs[diri];

				// Starts with 'processors'
				if (d.find("processors") == 0)
				{
					if (fileHandler().exists(d))
					{
						fileHandler().rmDir(d);
					}
				}

				// Starts with 'processor'
				if (d.find("processor") == 0)
				{
					// Check that integer after processor
					fileName num(d.substr(9));
					label proci = -1;
					if (Foam::read(num.c_str(), proci))
					{
						if (fileHandler().exists(d))
						{
							fileHandler().rmDir(d);
						}
					}
				}
			}

			procDirsProblem = false;
		}

		if (procDirsProblem)
		{
			FatalErrorInFunction
				<< "Case is already decomposed with " << nProcs
				<< " domains, use the -force option or manually" << nl
				<< "remove processor directories before decomposing. e.g.,"
				<< nl
				<< "    rm -rf " << runTime.path().c_str() << "/processor*"
				<< nl
				<< exit(FatalError);
		}
	}

    fvSolution solDict(runTime);
    dictionary harmonics = solDict.subDict("harmonicBalance");
    int nO = harmonics.getOrDefault<label>("instantsNumber",3);

    Info<< "Create mesh" << endl;
    domainDecomposition mesh
    (
        IOobject
        (
			fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        decompDictFile
    );

    // Decompose the mesh
    if (!decomposeFieldsOnly)
    {
        mesh.decomposeMesh();

        mesh.writeDecomposition(decomposeSets);

        if (writeCellDist)
        {
            const labelList& procIds = mesh.cellToProc();

            // Write decomposition as volScalarField for visualization
            volScalarField cellDist
            (
                IOobject
                (
                    "cellDist",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("cellDist", dimless, -1),
                zeroGradientFvPatchScalarField::typeName
            );

            forAll(procIds, celli)
            {
               cellDist[celli] = procIds[celli];
            }

            cellDist.correctBoundaryConditions();
            cellDist.write();

            Info<< nl << "Wrote decomposition as volScalarField to "
                << cellDist.name() << " for visualization."
                << endl;

            // Write decomposition as labelList for use with 'manual'
            // decomposition method.
            labelIOList cellDecomposition
            (
                IOobject
                (
                    "cellDecomposition",
                    mesh.facesInstance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procIds
            );
            cellDecomposition.write();

            Info<< nl << "Wrote decomposition to "
                << cellDecomposition.objectPath()
                << " for use in manual decomposition." << endl;
        }

        fileHandler().flush();

    	if(fileHandler().exists(cwd()/"tmpCnst"))
    	{
    		rmDir(cwd()/"tmpCnst");
    	}

    	if(fileHandler().exists(cwd()/"tmpSys"))
    	{
    		rmDir(cwd()/"tmpSys");
    	}

        for (int subLeveli = 0; subLeveli<nO; subLeveli++)
        {
        	fileName subLevConst = runTime.constant()/"subTimeLevel" + Foam::name(subLeveli);

        	if(fileHandler().exists(subLevConst))
        	{
        		rmDir(subLevConst);
        	}

    		fileName directorySys = runTime.system()/"subTimeLevel" + Foam::name(subLeveli);

    		if(fileHandler().exists(directorySys))
    		{
    			rmDir(directorySys);
    		}
        }

    	cp(runTime.system(), "tmpSys");
    	cp(runTime.constant(), "tmpCnst");

        for (int subLeveli = 0; subLeveli<nO; subLeveli++)
        {
    		cp("tmpCnst",runTime.constant()/"subTimeLevel" + Foam::name(subLeveli));
    		cp("tmpSys",runTime.system()/"subTimeLevel" + Foam::name(subLeveli));
        }

        fileName thermo = runTime.constant()/"thermophysicalProperties";
		fileName turbulence = runTime.constant()/"turbulenceProperties";

    	for (label proci = 0; proci < mesh.nProcs(); ++proci)
    	{
    		cp(cwd()/"processor" + Foam::name(proci)/"constant", cwd()/"processor" + Foam::name(proci)/"tmpCnst");

    		fileName procCnstDir = cwd()/"processor" + Foam::name(proci)/"tmpCnst";

    		for (int subLeveli = 0; subLeveli<nO; subLeveli++)
    		{
    			fileName procSubTimeDir =
    					cwd()/"processor" + Foam::name(proci)/"constant"/"subTimeLevel" + Foam::name(subLeveli);

    			cp(procCnstDir,procSubTimeDir);

    			cp(thermo, procSubTimeDir);
    			cp(turbulence, procSubTimeDir);
    		}

    		rmDir(procCnstDir);
    	}

    	if(fileHandler().exists(cwd()/"tmpCnst"))
    	{
    		rmDir(cwd()/"tmpCnst");
    	}

    	if(fileHandler().exists(cwd()/"tmpSys"))
    	{
    		rmDir(cwd()/"tmpSys");
    	}
    }

	for (int subLeveli = 0; subLeveli<nO; subLeveli++)
	{
        const word& subLevelName = "subTimeLevel" + Foam::name(subLeveli);
        const word& subLevelDir = subLevelName;

        mesh.polyMesh::rename(subLevelName);

        if (copyZero)
        {
            // Copy the 0 directory into each of the processor directories
            fileName prevTimePath;
            for (label proci = 0; proci < mesh.nProcs(); ++proci)
            {
                Time processorDb
                (
                    Time::controlDictName,
                    args.rootPath(),
                    args.caseName()/("processor" + Foam::name(proci))
                );
                processorDb.setTime(runTime);

                if (fileHandler().isDir(runTime.timePath()))
                {
                    // Get corresponding directory name (to handle processors/)
                    const fileName timePath
                    (
                        fileHandler().objectPath
                        (
                            IOobject
                            (
                                "",
                                processorDb.timeName(),
                                processorDb
                            ),
                            word::null
                        )
                    );

                    if (timePath != prevTimePath)
                    {
                        Info<< "Processor " << proci
                            << ": copying " << runTime.timePath() << nl
                            << " to " << timePath << endl;
                        fileHandler().cp(runTime.timePath(), timePath);

                        prevTimePath = timePath;
                    }
                }
            }
        }
        else
        {
            // Decompose the field files

            // Cached processor meshes and maps. These are only preserved if
            // running with multiple times.
            PtrList<Time> processorDbList(mesh.nProcs());
            PtrList<fvMesh> procMeshList(mesh.nProcs());
            PtrList<labelIOList> faceProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> cellProcAddressingList(mesh.nProcs());
            PtrList<labelIOList> boundaryProcAddressingList(mesh.nProcs());
            PtrList<fvFieldDecomposer> fieldDecomposerList(mesh.nProcs());
            PtrList<dimFieldDecomposer> dimFieldDecomposerList(mesh.nProcs());
            PtrList<labelIOList> pointProcAddressingList(mesh.nProcs());
            PtrList<pointFieldDecomposer> pointFieldDecomposerList
            (
                mesh.nProcs()
            );


            // Loop over all times
            forAll(times, timeI)
            {
                runTime.setTime(times[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl;

                // Search for list of objects for this time
                IOobjectList objects(mesh, runTime.timeName());

                // Construct the vol fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~
                PtrList<volScalarField> volScalarFields;
                readFields(mesh, objects, volScalarFields, false);
                PtrList<volVectorField> volVectorFields;
                readFields(mesh, objects, volVectorFields, false);
                PtrList<volSphericalTensorField> volSphericalTensorFields;
                readFields(mesh, objects, volSphericalTensorFields, false);
                PtrList<volSymmTensorField> volSymmTensorFields;
                readFields(mesh, objects, volSymmTensorFields, false);
                PtrList<volTensorField> volTensorFields;
                readFields(mesh, objects, volTensorFields, false);



                // Construct the dimensioned fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
                readFields(mesh, objects, dimScalarFields);
                PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
                readFields(mesh, objects, dimVectorFields);
                PtrList<DimensionedField<sphericalTensor, volMesh>>
                    dimSphericalTensorFields;
                readFields(mesh, objects, dimSphericalTensorFields);
                PtrList<DimensionedField<symmTensor, volMesh>>
                    dimSymmTensorFields;
                readFields(mesh, objects, dimSymmTensorFields);
                PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;
                readFields(mesh, objects, dimTensorFields);

                // Construct the surface fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                PtrList<surfaceScalarField> surfaceScalarFields;
                readFields(mesh, objects, surfaceScalarFields, false);
                PtrList<surfaceVectorField> surfaceVectorFields;
                readFields(mesh, objects, surfaceVectorFields, false);
                PtrList<surfaceSphericalTensorField>
                    surfaceSphericalTensorFields;
                readFields(mesh, objects, surfaceSphericalTensorFields, false);
                PtrList<surfaceSymmTensorField> surfaceSymmTensorFields;
                readFields(mesh, objects, surfaceSymmTensorFields, false);
                PtrList<surfaceTensorField> surfaceTensorFields;
                readFields(mesh, objects, surfaceTensorFields, false);

                // Construct the point fields
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~
                const pointMesh& pMesh = pointMesh::New(mesh);

                PtrList<pointScalarField> pointScalarFields;
                readFields(pMesh, objects, pointScalarFields, false);
                PtrList<pointVectorField> pointVectorFields;
                readFields(pMesh, objects, pointVectorFields, false);
                PtrList<pointSphericalTensorField> pointSphericalTensorFields;
                readFields(pMesh, objects, pointSphericalTensorFields, false);
                PtrList<pointSymmTensorField> pointSymmTensorFields;
                readFields(pMesh, objects, pointSymmTensorFields, false);
                PtrList<pointTensorField> pointTensorFields;
                readFields(pMesh, objects, pointTensorFields, false);

                Info<< endl;

                // split the fields over processors
                for (label proci = 0; proci < mesh.nProcs(); ++proci)
                {
                    Info<< "Processor " << proci << ": field transfer" << endl;

                    // open the database
                    if (!processorDbList.set(proci))
                    {
                        processorDbList.set
                        (
                            proci,
                            new Time
                            (
                                Time::controlDictName,
                                args.rootPath(),
                                args.caseName()
                              / ("processor" + Foam::name(proci))
                            )
                        );
                    }
                    Time& processorDb = processorDbList[proci];

                    processorDb.setTime(runTime);

                    // read the mesh
                    if (!procMeshList.set(proci))
                    {
                        procMeshList.set
                        (
                            proci,
                            new fvMesh
                            (
                                IOobject
                                (
									subLevelName,
                                    processorDb.timeName(),
                                    processorDb
                                )
                            )
                        );
                    }

                    const fvMesh& procMesh = procMeshList[proci];

                    const labelIOList& faceProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "faceProcAddressing",
                        faceProcAddressingList
                    );

                    const labelIOList& cellProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "cellProcAddressing",
                        cellProcAddressingList
                    );

                    const labelIOList& boundaryProcAddressing = procAddressing
                    (
                        procMeshList,
                        proci,
                        "boundaryProcAddressing",
                        boundaryProcAddressingList
                    );

                    // FV fields
                    {
                        if (!fieldDecomposerList.set(proci))
                        {
                            fieldDecomposerList.set
                            (
                                proci,
                                new fvFieldDecomposer
                                (
                                    mesh,
                                    procMesh,
                                    faceProcAddressing,
                                    cellProcAddressing,
                                    boundaryProcAddressing
                                )
                            );
                        }
                        const fvFieldDecomposer& fieldDecomposer =
                            fieldDecomposerList[proci];

                        fieldDecomposer.decomposeFields(volScalarFields);
                        fieldDecomposer.decomposeFields(volVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            volSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields(volSymmTensorFields);
                        fieldDecomposer.decomposeFields(volTensorFields);

                        fieldDecomposer.decomposeFields(surfaceScalarFields);
                        fieldDecomposer.decomposeFields(surfaceVectorFields);
                        fieldDecomposer.decomposeFields
                        (
                            surfaceSphericalTensorFields
                        );
                        fieldDecomposer.decomposeFields
                        (
                            surfaceSymmTensorFields
                        );
                        fieldDecomposer.decomposeFields(surfaceTensorFields);

                        if (times.size() == 1)
                        {
                            // Clear cached decomposer
                            fieldDecomposerList.set(proci, nullptr);
                        }
                    }

                    // Dimensioned fields
                    {
                        if (!dimFieldDecomposerList.set(proci))
                        {
                            dimFieldDecomposerList.set
                            (
                                proci,
                                new dimFieldDecomposer
                                (
                                    mesh,
                                    procMesh,
                                    faceProcAddressing,
                                    cellProcAddressing
                                )
                            );
                        }
                        const dimFieldDecomposer& dimDecomposer =
                            dimFieldDecomposerList[proci];

                        dimDecomposer.decomposeFields(dimScalarFields);
                        dimDecomposer.decomposeFields(dimVectorFields);
                        dimDecomposer.decomposeFields(dimSphericalTensorFields);
                        dimDecomposer.decomposeFields(dimSymmTensorFields);
                        dimDecomposer.decomposeFields(dimTensorFields);

                        if (times.size() == 1)
                        {
                            dimFieldDecomposerList.set(proci, nullptr);
                        }
                    }

                    // Point fields
                    if
                    (
                        pointScalarFields.size()
                     || pointVectorFields.size()
                     || pointSphericalTensorFields.size()
                     || pointSymmTensorFields.size()
                     || pointTensorFields.size()
                    )
                    {
                        const labelIOList& pointProcAddressing = procAddressing
                        (
                            procMeshList,
                            proci,
                            "pointProcAddressing",
                            pointProcAddressingList
                        );

                        const pointMesh& procPMesh = pointMesh::New(procMesh);

                        if (!pointFieldDecomposerList.set(proci))
                        {
                            pointFieldDecomposerList.set
                            (
                                proci,
                                new pointFieldDecomposer
                                (
                                    pMesh,
                                    procPMesh,
                                    pointProcAddressing,
                                    boundaryProcAddressing
                                )
                            );
                        }
                        const pointFieldDecomposer& pointDecomposer =
                            pointFieldDecomposerList[proci];

                        pointDecomposer.decomposeFields(pointScalarFields);
                        pointDecomposer.decomposeFields(pointVectorFields);
                        pointDecomposer.decomposeFields
                        (
                            pointSphericalTensorFields
                        );
                        pointDecomposer.decomposeFields(pointSymmTensorFields);
                        pointDecomposer.decomposeFields(pointTensorFields);


                        if (times.size() == 1)
                        {
                            pointProcAddressingList.set(proci, nullptr);
                            pointFieldDecomposerList.set(proci, nullptr);
                        }
                    }

                    // Decompose the "uniform" directory in the time region
                    // directory
                    decomposeUniform(copyUniform, mesh, processorDb, subLevelDir);

					if (subLeveli == 0)
					{
						decomposeUniform(copyUniform, mesh, processorDb);
					}

                    // We have cached all the constant mesh data for the current
                    // processor. This is only important if running with
                    // multiple times, otherwise it is just extra storage.
                    if (times.size() == 1)
                    {
                        boundaryProcAddressingList.set(proci, nullptr);
                        cellProcAddressingList.set(proci, nullptr);
                        faceProcAddressingList.set(proci, nullptr);
                        procMeshList.set(proci, nullptr);
                        processorDbList.set(proci, nullptr);
                    }
                }
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
