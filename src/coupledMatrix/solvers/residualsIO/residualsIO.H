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

Class
    Foam::residualIO

Description
    Small wrapper class to input and output residual norms.

Authors

\*---------------------------------------------------------------------------*/

#ifndef residualsIO_H
#define residualsIO_H

#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class residualIO Declaration
\*---------------------------------------------------------------------------*/

class residualsIO
{
	// Private data

		//- Initial residuals
		List<scalar> sInitRes_;
		List<vector> vInitRes_;

		//- Final residuals
		List<scalar> sFinalRes_;
		List<vector> vFinalRes_;

		//- Number of iterations
		label nIterations_;

		//- Has the solver converged?
		bool converged_;

public:

    residualsIO(const label nScalars, const label nVectors)
     :
        sInitRes_(nScalars, 0.0),
        vInitRes_(nVectors, Zero),
        sFinalRes_(nScalars, 0.0),
        vFinalRes_(nVectors, Zero),
 		nIterations_(0),
 		converged_(false)
    {}

    //- Initialise with uniform values
    residualsIO(const residualsIO& res)
     :
        sInitRes_(res.sInitRes_),
        vInitRes_(res.vInitRes_),
        sFinalRes_(res.sFinalRes_),
        vFinalRes_(res.vFinalRes_),
		nIterations_(res.nIterations_),
		converged_(res.converged_)
    {}

	//- Return initial residual
	scalarList& sInitRes()
	{
		return sInitRes_;
	}

	List<vector>& vInitRes()
	{
		return vInitRes_;
	}

	//- Return final residual
	scalarList& sFinalRes()
	{
		return sFinalRes_;
	}

	List<vector>& vFinalRes()
	{
		return vFinalRes_;
	}

	//- Return number of iterations
	label& nIterations()
	{
		return nIterations_;
	}

	//- Has the solver converged?
	bool& converged()
	{
		return converged_;
	}

	//- Print summary of solver performance
	void print() const
	{
		Info<<"Initial residual coupled solver = ";

		forAll(sInitRes_, j)
		{
			Info << " ";
			Info << sInitRes_[j];
		}
		forAll(vInitRes_, j)
		{
			Info << " ";
			Info << vInitRes_[j];
		}

		Info<<endl;

		Info<<"Final residual coupled solver = ";

		forAll(sFinalRes_, j)
		{
			Info << " ";
			Info << sFinalRes_[j];
		}
		forAll(vFinalRes_, j)
		{
			Info << " ";
			Info << vFinalRes_[j];
		}

		Info<<endl;

		Info << "No Iterations coupled solver " << nIterations_
		<< endl;
	}


	const residualsIO& operator=(const residualsIO& res)
    {
        sInitRes_  = res.sInitRes_;
        vInitRes_ = res.vInitRes_;

        sFinalRes_  = res.sFinalRes_;
        vFinalRes_ = res.vFinalRes_;

        nIterations_ = res.nIterations_;

        converged_ = res.converged_;

        return *this;
    }


    friend residualsIO operator/(const residualsIO& val1, const residualsIO& val2)
    {
        residualsIO res(val1);

        forAll(val2.sInitRes_, i)
        {
            res.sInitRes_[i] /= val2.sInitRes_[i] + ROOTVSMALL;
            res.sFinalRes_[i] /= val2.sFinalRes_[i] + ROOTVSMALL;
        }

        vector smallVector(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        forAll(val2.vInitRes_, i)
        {
            res.vInitRes_[i] = cmptDivide(res.vInitRes_[i], val2.vInitRes_[i] + smallVector);
            res.vFinalRes_[i] = cmptDivide(res.vFinalRes_[i], val2.vFinalRes_[i] + smallVector);
        }

        return res;
    }

	scalar maxRel(vector::labelType validComponents) const
	{
        scalar maxValRel(-VGREAT);

        forAll(sFinalRes_, i)
        {
            maxValRel = max(maxValRel,sFinalRes_[i]/(sInitRes_[i] + ROOTVSMALL));
        }
        forAll(vFinalRes_, i)
        {
            forAll(validComponents, cmpt)
            {
                if (validComponents[cmpt] == 1)
                {
                	scalar cmptRelRes = vFinalRes_[i][cmpt]/(vInitRes_[i][cmpt] + ROOTVSMALL);
                	maxValRel = max(maxValRel, cmptRelRes);
                }
            }
        }

        return maxValRel;
	}

    friend scalar max(const residualsIO& res)
    {
        scalar maxVal(-VGREAT);

        forAll(res.sFinalRes_, i)
        {
            maxVal = max(maxVal,res.sFinalRes_[i]);
        }
        forAll(res.vFinalRes_, i)
        {
            maxVal = max(maxVal,cmptMax(res.vFinalRes_[i]));
        }
        return maxVal;
    }


    friend scalar maxInit(const residualsIO& res)
    {
        scalar maxVal(-VGREAT);

        forAll(res.sInitRes_, i)
        {
            maxVal = max(maxVal,res.sInitRes_[i]);
        }
        forAll(res.vInitRes_, i)
        {
            maxVal = max(maxVal,cmptMax(res.vInitRes_[i]));
        }
        return maxVal;
    }

    friend Istream& operator>>(Istream& is, residualsIO& res)
    {

		forAll(res.sInitRes_, i)
		{
			is >> res.sInitRes_[i];
		}

		forAll(res.vInitRes_, i)
		{
			is >> res.vInitRes_[i];
		}

		forAll(res.sFinalRes_, i)
		{
			is >> res.sFinalRes_[i];
		}

		forAll(res.vFinalRes_, i)
		{
			is >> res.vFinalRes_[i];
		}

        return is;
    }

    friend Ostream& operator<<(Ostream& os, const residualsIO& res)
    {
		forAll(res.sInitRes_, i)
		{
			os << res.sInitRes_[i];
			os << " ";
		}

		forAll(res.vInitRes_, i)
		{
			os << res.vInitRes_[i];
			os << " ";
		}

		forAll(res.sFinalRes_, i)
		{
			os << res.sFinalRes_[i];
			os << " ";
		}

		forAll(res.vFinalRes_, i)
		{
			os << res.vFinalRes_[i];
			os << " ";
		}

        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
