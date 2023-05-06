/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2015 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2015 Oliver Oxtoby - CSIR, South Africa

    Copyright (C) 2022 Stefano Oliani
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

\*---------------------------------------------------------------------------*/

#include "pseudotimeControl.H"
#include "Switch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudotimeControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


bool Foam::pseudotimeControl::read()
{
    solutionControl::read(false);

    // Read solution controls
    const dictionary& pseudoDict = dict();

    nCorrOuter_ = pseudoDict.lookupOrDefault<label>("nPseudoCorr", 20);
    nCorrOuterMin_ = pseudoDict.lookupOrDefault<label>("nPseudoCorrMin", 1);

    pseudoDict.lookup("pseudoTol") >> residualTols_;

    pseudoDict.lookup("pseudoTolRel") >> residualTolsRel_;

    turbOnFinalIterOnly_ =
        pseudoDict.getOrDefault("turbOnFinalIterOnly", false);

    if (state_.found("initResiduals"))
    {
        state_.lookup("initResiduals") >> initResiduals_;
    }
    else
    {
        initResiduals_ = residualsIO(residuals_.sInitRes().size(), residuals_.vInitRes().size());
    }

    return true;
}


bool Foam::pseudotimeControl::criteriaSatisfied()
{
    // no checks on first iteration after restart or first iteration of 
    // timestep - nothing has been calculated yet
    if (firstIteration_ || corr_ == 1 || finalIter())
    {
        firstIteration_ = false;
        return false;
    }

    bool storeIni = this->storeInitialResiduals();

    if (storeIni)
    {
        initResiduals_ = residuals_;

        state_.add("initResiduals", initResiduals_, true); // Store for restarts
    }

    bool absCheck = (maxInit(residuals_) < residualTols_);

    bool relCheck = false;

    if (!storeIni)
    {
        relCheck = (maxInit(residuals_/(initResiduals_)) < residualTolsRel_);
    }

    return (corr_ >= nCorrOuterMin_) && (absCheck || relCheck);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pseudotimeControl::pseudotimeControl
(
	fvMesh& mesh,
	const bool steadyState,
	const label nScalars,
	const label nVectors
)
:
    solutionControl(mesh, "pseudoTime"),
    steadyState_(steadyState),
    nCorrOuter_(0),
    turbOnFinalIterOnly_(false),
    converged_(false),
    firstIteration_(true),
    residualTols_(),
    residualTolsRel_(),
    residuals_(nScalars, nVectors),
	initResiduals_(nScalars, nVectors),
    state_
    (
        IOobject
        (
            "pseudotimeState",
            mesh.time().timeName(),
            "uniform",
            mesh.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    )
{
    read();

    Info<< algorithmName_ << ": ";
    if (!steadyState)
    {
        Info<< "max iterations = " << nCorrOuter_ << ", ";
    }
    Info<< "tolerance = " << residualTols_ << ", "
        << "relTol = " << residualTolsRel_
        << nl;
    Info<< endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudotimeControl::~pseudotimeControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pseudotimeControl::loop()
{
    read();

    corr_++;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corr = " << corr_ << endl;
    }

    // If steady state, ignore nCorrOuter because endTime is used 
    // to terminate main loop
    if (!steadyState_ && corr_ == nCorrOuter_ + 1)
    {
        if (nCorrOuter_ != 1)
        {
            Info<< algorithmName_ << ": not converged within "
                << nCorrOuter_ << " iterations" << endl;
        }

        corr_ = 0;
        if (!steadyState_)
        {
            mesh_.data::remove("finalIteration");
        }
        return false;
    }

    bool completed = false;
    if (converged_ || criteriaSatisfied())
    {
        if (converged_)
        {
            Info<< algorithmName_ << ": converged in " << corr_ - 1
                << " iterations" << endl;

            if (!steadyState_)
            {
                mesh_.data::remove("finalIteration");
            }
            corr_ = 0;
            converged_ = false;

            completed = true;
        }
        else
        {
            Info<< algorithmName_ << ": iteration " << corr_ << endl;
            storePrevIterFields();

            if (!steadyState_)
            {
                mesh_.data::add("finalIteration", true);
            }
            converged_ = true;
        }
    }
    else
    {
        if (finalIter())
        {
            if (!steadyState_)
            {
                mesh_.data::add("finalIteration", true);
            }
        }

        if (steadyState_ || corr_ <= nCorrOuter_)
        {
            if (nCorrOuter_ != 1)
            {
                Info<< algorithmName_ << ": iteration " << corr_ << endl;

                storePrevIterFields();
            }

            completed = false;
        }
    }

    return !completed;
}


// ************************************************************************* //
