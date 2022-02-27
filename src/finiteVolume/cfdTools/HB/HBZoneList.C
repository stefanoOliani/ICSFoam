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

#include "HBZoneList.H"

#include "volFields.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::HBZoneList::findT0() const
{
	scalar minOmega = GREAT;
	scalar minOmegai = GREAT;

    forAll(*this, i)
    {
        const scalarList& omegaListi = operator[](i).omegaList();

        forAll (omegaListi, k)
        {
        	const scalar& omegaAbsVal = std::abs(omegaListi[k]);

        	if (omegaAbsVal > 0 && omegaAbsVal < minOmegai)
        	{
        		minOmegai = omegaAbsVal;
        	}
        }

        if (minOmegai < minOmega)
        {
        	minOmega = minOmegai;
        }
    }

    return 2*constant::mathematical::pi/minOmega;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HBZoneList::HBZoneList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<HBZone>(),
	selectedSnapshots_(1,0.0),
    mesh_(mesh)
{
    reset(dict);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::HBZoneList::~HBZoneList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::HBZoneList::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "    No HB zones active" << endl;
    }

    return a;
}


void Foam::HBZoneList::reset(const dictionary& dict)
{
    label count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            ++count;
        }
    }

    this->resize(count);

    count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            const word& name = dEntry.keyword();
            const dictionary& modelDict = dEntry.dict();

            Info<< "    creating HB zone: " << name << endl;

            this->set
            (
                count++,
                new HBZone(name, mesh_, modelDict)
            );
        }
    }
}


const Foam::scalarList& Foam::HBZoneList::selectedSnapshots() const
{
	return selectedSnapshots_;
}


bool Foam::HBZoneList::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        HBZone& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::HBZoneList::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::HBZoneList::addingSource
(
	coupledMatrix& eqSystem,
	PtrList<PtrList<volScalarField>>& scalarVars,
	PtrList<PtrList<volVectorField>>& vectorVars,
	label instantNo
) const
{
    forAll(*this, i)
    {
        operator[](i).addSource(eqSystem, scalarVars, vectorVars, instantNo);
    }
}


void Foam::HBZoneList::addingSource
(
	PtrList<vectorField>& source,
	PtrList<PtrList<volVectorField>>& vars,
	label varNum,
	bool cylCoords
) const
{
    forAll(*this, i)
    {
        operator[](i).addSource(source, vars, varNum, cylCoords);
    }
}


void Foam::HBZoneList::addSource
(
	fvVectorMatrix& eqn,
	PtrList<volVectorField>& fieldPtr,
	label instantNo,
	bool  cylCoords
) const
{
    forAll(*this, i)
    {
        operator[](i).addSource(eqn, fieldPtr, instantNo, cylCoords);
    }
}

void Foam::HBZoneList::addBlock
(
	blockFvMatrix<scalar,scalar>& blockMatrix,
	label rowIndex,
	label colIndex
) const
{
    forAll(*this, i)
    {
        operator[](i).addBlock(blockMatrix, rowIndex, colIndex);
    }
}


void Foam::HBZoneList::addBlock
(
	blockFvMatrix<vector,tensor>& blockMatrix,
	label rowIndex,
	label colIndex
) const
{
    forAll(*this, i)
    {
        operator[](i).addBlock(blockMatrix, rowIndex, colIndex);
    }
}


void Foam::HBZoneList::setInstants()
{
	const scalar T0 = findT0();

    const dictionary& solDict = mesh_.solutionDict();

	const label nInstants =
			solDict.subDict("harmonicBalance").getOrDefault<label>("instantsNumber",3);

	if (solDict.subDict("harmonicBalance").found("selectedPeriod"))
	{
		scalar selT = solDict.subDict("harmonicBalance").get<scalar>("selectedPeriod");
		scalar dT = selT/nInstants;

		scalarList selSnapshots(nInstants, 0.0);

		forAll (selSnapshots, snapi)
		{
			selSnapshots[snapi] = snapi*dT;
		}

		selectedSnapshots_ = selSnapshots;

		Info<<"The selected snapshots are "<<
				selSnapshots<<endl;
	}
	else
	{
		scalar Tf = 5*T0;
		const scalar stepSize = 0.001*T0;

		scalarList timePeriodsSet(1,T0);
		scalar Tfi = T0;

		while (Tfi <= Tf)
		{
			Tfi += stepSize;
			timePeriodsSet.append(Tfi);
		}

		scalarList snapshots(nInstants, 0.0);
		scalarList selSnapshots(nInstants, 0.0);
		scalar condNumber = GREAT;

		forAll(timePeriodsSet, k)
		{
			scalar TFinal = timePeriodsSet[k];
			scalar deltaT = TFinal/nInstants;

			forAll(snapshots, l)
			{
				snapshots[l] = deltaT*l;
			}

			scalar maxZoneCondN = SMALL;

			forAll(*this, i)
			{
				scalar zoneCondNi = operator[](i).calcConditionN(snapshots);

				if (zoneCondNi > maxZoneCondN)
				{
					maxZoneCondN = zoneCondNi;
				}
			}

			if (maxZoneCondN < condNumber)
			{
				condNumber = maxZoneCondN;
				selSnapshots = snapshots;
			}
		}

		selectedSnapshots_ = selSnapshots;

		Info<<"The maximum condition number among the HBZones is "<<
				condNumber<<endl;

		if (condNumber > 4.0)
		{
			WarningInFunction<<
					"The condition number in the search range is high"<<endl;
		}

		Info<<"The selected snapshots are "<<
				selSnapshots<<endl;
	}

    forAll(*this, i)
    {
        operator[](i).updateHBOperators(selectedSnapshots_);
    }

}


void Foam::HBZoneList::factStep
(
	PtrList<PtrList<volScalarField>>& scalarVars,
	PtrList<PtrList<volVectorField>>& vectorVars,
	PtrList<volScalarField>& deltaTField
) const
{
    forAll(*this, i)
    {
        operator[](i).factorizationStep(scalarVars, vectorVars, deltaTField);
    }
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const HBZoneList& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
