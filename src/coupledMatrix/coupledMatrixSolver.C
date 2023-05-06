/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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

#include "error.H"

#include "coupledMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

	defineRunTimeSelectionTable(coupledMatrix::solver, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<coupledMatrix::solver> coupledMatrix::solver::New
(
    const dictionary& dict,
    const coupledMatrix& matrix
)
{
    word solverTypeName;

    dict.lookup("solver") >> solverTypeName;

    auto* ctorPtr = dictionaryConstructorTable(solverTypeName);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown solver type " << solverTypeName
            << endl << endl
            << "Valid solver types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solver>
           (
               ctorPtr(dict.subDict(solverTypeName), matrix)
           );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledMatrix::solver::solver
(
     const word& type,
     const dictionary& dict,
     const coupledMatrix& matrix
)
:
     matrix_(matrix),
     controlDict_(dict),
	 maxIter_(controlDict_.getOrDefault<label>("maxIter", defaultMaxIter_)),
	 minIter_(controlDict_.getOrDefault<label>("minIter", 0)),
	 tolerance_(controlDict_.get<scalar>("tolerance")),
	 relTolerance_(controlDict_.get<scalar>("relTol"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMatrix::solver::readControls()
{
    maxIter_ = controlDict_.getOrDefault<label>("maxIter", defaultMaxIter_);
    minIter_ = controlDict_.getOrDefault<label>("minIter", 0);
    tolerance_ = controlDict_.get<scalar>("tolerance");
    relTolerance_ = controlDict_.get<scalar>("relTol");
}


void Foam::coupledMatrix::solver::normFactors
(
	PtrList<volScalarField>& sW,
	PtrList<volVectorField>& vW,
	const PtrList<scalarField>& sSource,
	const PtrList<vectorField>& vSource,
	scalarList& sNormFactor,
	List<vector>& vNormFactor
) const
{
    const coupledMatrix& cMatrix = this->matrix();
    const int nScalar = cMatrix.nScal();
    const int nVector = cMatrix.nVect();

    PtrList<scalarField> wsTmp(nScalar);
    PtrList<vectorField> wvTmp(nVector);
    PtrList<scalarField> psTmp(nScalar);
    PtrList<vectorField> pvTmp(nVector);

    scalarList sRef(nScalar);
    List<vector> vRef(nVector);

    PtrList<volScalarField> sRefField(nScalar);
    PtrList<volVectorField> vRefField(nVector);

    forN(nScalar,j)
    {
        wsTmp.set(j, new scalarField(sW[j].size(), 0.0));
        psTmp.set(j, new scalarField(sW[j].size(), 0.0));

        sRef[j] = gAverage(sW[j]);

        sRefField.set
        (
            j,
            new volScalarField
            (
                IOobject
                (
                    "sRefField" + sW[j].name(),
                    matrix_.mesh().time().timeName(),
                    matrix_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				matrix_.mesh(),
				dimensionedScalar("sRefField", dimless, sRef[j])
            )
        );

    }
    forN(nVector,j)
    {
    	wvTmp.set(j, new vectorField(vW[j].size(), Zero));
    	pvTmp.set(j, new vectorField(vW[j].size(), Zero));

    	vRef[j] = gAverage(vW[j]);


        vRefField.set
        (
            j,
            new volVectorField
            (
                IOobject
                (
                    "vRefField" + vW[j].name(),
                    matrix_.mesh().time().timeName(),
                    matrix_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
				matrix_.mesh(),
				dimensionedVector("vRefField", dimless, vRef[j])
            )
        );
    }

    // Calculate A.x
    this->matrix_.matrixMul(sW, vW, wsTmp, wvTmp);

    // Calculate A.xRef
    this->matrix_.matrixMul(sRefField, vRefField, psTmp, pvTmp);

    forAll(sNormFactor, i)
    {
        sNormFactor[i] = gSum(mag(wsTmp[i] - psTmp[i]) + mag(sSource[i]- psTmp[i])) + SMALL;
    }

    vector small = pTraits<vector>::one*SMALL;

    forAll(vNormFactor, i)
    {
        vNormFactor[i] = gSum(cmptMag(wvTmp[i] - pvTmp[i]) + cmptMag(vSource[i] - pvTmp[i])) + small;
    }

}


bool Foam::coupledMatrix::solver::stop
(
    residualsIO& solverPerf
) const
{
    if (solverPerf.nIterations() < minIter_)
    {
        return false;
    }

    // Make sure the relative tolerance is ignored in non-solved directions
    vector::labelType validComponents = this->matrix().mesh().solutionD(); //-1 for empty directions

    if (   (solverPerf.nIterations() >= maxIter_)
		|| (max(solverPerf) < tolerance_)
		|| (solverPerf.maxRel(validComponents) < relTolerance_)
	   )
    {
        return true;
    }

    return false;

}

void Foam::coupledMatrix::solver::read(const dictionary& solverControls)
{
    controlDict_ = solverControls;
    readControls();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
