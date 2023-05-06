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
    jacobianMatrix

Description
    Holds the different components of a coupled matrix in one place.

Authors


\*---------------------------------------------------------------------------*/

#ifndef coupledMatrix_H
#define coupledMatrix_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "blockFvMatrix.H"
#include "blockFvOperators.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "residualsIO.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forN(N,i) for(Foam::label i=0; i < (N); i++)

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class coupledMatrix Declaration
\*---------------------------------------------------------------------------*/

class coupledMatrix
{
		//- Number of scalar and vector unknowns
		int	nScal_, nVect_;

		//- Is the system in Delta form?
		bool deltaForm_;


    protected:

        // Components of coupled matrix

		//- Effect of scalars on scalars
        PtrList<blockFvMatrix<scalar,scalar>> dSByS_;

        //- Effect of vectors on scalars
        PtrList<blockFvMatrix<scalar,vector>> dSByV_;

        //- Effect of scalars on vectors
        PtrList<blockFvMatrix<vector,vector>> dVByS_;

        //- Effect of vectors on vectors
        PtrList<blockFvMatrix<vector,tensor>> dVByV_;

        const fvMesh& mesh_;

    public:

        //- Abstract base-class for LduMatrix solvers
            class solver
            {
            protected:

                // Protected data

                    const coupledMatrix& matrix_;

                    //- Dictionary of controls
                    dictionary controlDict_;

                    //- Default maximum number of iterations in the solver
                    static const label defaultMaxIter_;

                    //- Maximum number of iterations in the solver
                    label maxIter_;

                    //- Minimum number of iterations in the solver
                    label minIter_;

                    scalar tolerance_;

                    scalar relTolerance_;

                // Protected Member Functions

                    //- Read the control parameters from the controlDict_
                    virtual void readControls();

                    //- Compute normalisation factors
                    void normFactors
                    (
						PtrList<volScalarField>& sW,
						PtrList<volVectorField>& vW,
						const PtrList<scalarField>& sSource,
						const PtrList<vectorField>& vSource,
						scalarList& sNormFactor,
						List<vector>& vNormFactor
                    ) const;

            public:

                //- Runtime type information
                virtual const word& type() const = 0;

                // Declare run-time constructor selection table

				declareRunTimeSelectionTable
				(
					autoPtr,
					solver,
					dictionary,
					(
						const dictionary& dict,
						const coupledMatrix& matrix
					),
					(dict, matrix)
				);

                // Constructors

                solver
                (
                    const word& type,
                    const dictionary& dict,
                    const coupledMatrix& matrix
                );


                // Selectors

                    //- Return a new solver
                    static autoPtr<solver> New
                    (
                    	const dictionary& solverDict,
                        const coupledMatrix& matrix
                    );


                // Destructor

                    virtual ~solver() = default;


                // Member functions

                    // Access

                        const coupledMatrix& matrix() const
                        {
                            return matrix_;
                        }


                    //- Read and reset the solver parameters from the given dictionary
                    virtual void read(const dictionary& solverDict);

                    bool stop(residualsIO& solverPerf) const;

                    //- Solve
                    virtual residualsIO solve
                    (
                        PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
                        const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
                    ) const = 0;

                    //- Solve a system in Delta form
                    virtual residualsIO solveDelta
                    (
                        PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
                        const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource
                    ) const = 0;

                    virtual residualsIO solveDelta
                    (
                        PtrList<volScalarField>& sW, PtrList<volVectorField>& vW,
                        const PtrList<scalarField>& sSource, const PtrList<vectorField>& vSource,
            			PtrList<volScalarField>& dsW, PtrList<volVectorField>& dvW
                    ) const = 0;

            };


            //- Abstract base-class for LduMatrix smoothers
                class smoother
                {
                protected:

                    // Protected data

                        const coupledMatrix& matrix_;


                public:

                    //- Runtime type information
                    virtual const word& type() const = 0;


                    // Declare run-time constructor selection tables

                        declareRunTimeSelectionTable
                        (
                            autoPtr,
                            smoother,
                            dictionary,
                            (
								const coupledMatrix& matrix
                            ),
                            (
                                matrix
                            )
                        );


                    // Constructors

                        smoother
                        (
                            const coupledMatrix& matrix
                        );


                    // Selectors

                        //- Return a new smoother
                        static autoPtr<smoother> New
                        (
                            const coupledMatrix& matrix,
                            const dictionary& smootherDict
                        );


                    // Destructor

                        virtual ~smoother() = default;


                    // Member functions

                        // Access

                            const coupledMatrix& matrix() const
                            {
                                return matrix_;
                            }


                        //- Smooth the solution for a given number of sweeps
                        virtual void smooth
                        (
							PtrList<volScalarField>& sW,
							PtrList<volVectorField>& vW,
							const PtrList<scalarField>& sSource,
							const PtrList<vectorField>& vSource,
                            const label nSweeps
                        ) const = 0;
                };



            //- Abstract base-class for lduMatrix preconditioners
                class preconditioner
                {
                protected:

                    // Protected data

                        //- Reference to the base-solver this preconditioner is used with
                        const solver& solver_;


                public:

                    //- Find the preconditioner name (directly or from a sub-dictionary)
                    static word getName(const dictionary&);

                    //- Runtime type information
                    virtual const word& type() const = 0;


                    // Declare run-time constructor selection table

                        declareRunTimeSelectionTable
                        (
                            autoPtr,
                            preconditioner,
                            dictionary,
                            (
                                const solver& sol,
                                const dictionary& dict
                            ),
                            (sol, dict)
                        );


                    // Constructors

                        preconditioner
                        (
                            const solver& sol
                        )
                        :
                            solver_(sol)
                        {}


                    // Selectors

                        //- Return a pointer to the selected preconditioner
                        static autoPtr<preconditioner> New
                        (
                            const solver& sol,
                            const dictionary& dict
                        );


                    //- Destructor
                    virtual ~preconditioner() = default;


                    // Member functions

                        //- Read and reset the preconditioner parameters
                        //- from the given stream
                        virtual void read(const dictionary&)
                        {}

                        //- Precondition a vector
                        virtual void precondition
                        (
                            PtrList<scalarField>& sVec,
                            PtrList<vectorField>& vVec
                        ) const = 0;

                };




        //- Runtime type information
        TypeName("coupledMatrix");

        coupledMatrix
		(
			const fvMesh& mesh,
			int nScalars,
			int nVectors,
			bool deltaForm
		);

        virtual ~coupledMatrix()
        {}

        int nScal() const
        {
        	return nScal_;
        }

        int nVect() const
        {
        	return nVect_;
        }


        bool dSBySExists(const label i, const label j) const
        {
            const label idx = i*nScal_+j;
            return dSByS_.set(idx);
        }

        bool dSByVExists(const label i, const label j) const
        {
            const label idx = i*nVect_+j;
            return dSByV_.set(idx);
        }

        bool dVBySExists(const label i, const label j) const
        {
            const label idx = i*nScal_+j;
            return dVByS_.set(idx);
        }

        bool dVByVExists(const label i, const label j) const
        {
            const label idx = i*nVect_+j;
            return dVByV_.set(idx);
        }




        const blockFvMatrix<scalar,scalar>& dSByS(const label i, const label j) const
        {
            const label idx = i*nScal_+j;
            if (!dSByS_.set(idx))
            {
                FatalErrorInFunction
                    << "Coupled matrix component not defined." << exit(FatalError);
            }
            return dSByS_[idx];
        }

        const blockFvMatrix<scalar,vector>& dSByV(const label i, const label j) const
        {
            const label idx = i*nVect_+j;
            if (!dSByV_.set(idx))
            {
                FatalErrorInFunction
                    << "Coupled matrix component not defined." << exit(FatalError);
            }
            return dSByV_[idx];
        }

        const blockFvMatrix<vector,vector>& dVByS(const label i, const label j) const
        {
            const label idx = i*nScal_+j;
            if (!dVByS_.set(idx))
            {
                FatalErrorInFunction
                    << "Coupled matrix component not defined." << exit(FatalError);
            }
            return dVByS_[idx];
        }

        const blockFvMatrix<vector,tensor>& dVByV(const label i, const label j) const
        {
            const label idx = i*nVect_+j;
            if (!dVByV_.set(idx))
            {
                FatalErrorInFunction
                    << "Coupled matrix component not defined." << exit(FatalError);
            }
            return dVByV_[idx];
        }





        blockFvMatrix<scalar,scalar>& dSByS(const label i, const label j)
        {
            const label idx = i*nScal_+j;
            if (!dSByS_.set(idx))
            {
                dSByS_.set
                (
                    idx,
                    new blockFvMatrix<scalar,scalar>(mesh_)
                );
            }
            return dSByS_[idx];
        }

        blockFvMatrix<scalar,vector>& dSByV(const label i, const label j)
        {
            const label idx = i*nVect_+j;
            if (!dSByV_.set(idx))
            {
                dSByV_.set
                (
                    idx,
                    new blockFvMatrix<scalar,vector>(mesh_)
                );
            }
            return dSByV_[idx];
        }

        blockFvMatrix<vector,vector>& dVByS(const label i, const label j)
        {
            const label idx = i*nScal_+j;
            if (!dVByS_.set(idx))
            {
                dVByS_.set
                (
                    idx,
                    new blockFvMatrix<vector,vector>(mesh_)
                );
            }
            return dVByS_[idx];
        }

        blockFvMatrix<vector,tensor>& dVByV(const label i, const label j)
        {
            const label idx = i*nVect_+j;
            if (!dVByV_.set(idx))
            {
                dVByV_.set
                (
                    idx,
                    new blockFvMatrix<vector,tensor>(mesh_)
                );
            }
            return dVByV_[idx];
        }


        const fvMesh& mesh() const
        {
            return mesh_;
        }

        void matrixMul
        (
            PtrList<volScalarField>& sVec, PtrList<volVectorField>& vVec,
            PtrList<scalarField>& sResult, PtrList<vectorField>& vResult
        ) const;

        // Used in block Jacobi smoother
        void matrixMulNoDiag
		(
			PtrList<volScalarField>& sVec, PtrList<volVectorField>& vVec,
			PtrList<scalarField>& sResult, PtrList<vectorField>& vResult
		) const;

        // Solver calls for coupledMatrix

            //- Solve returning the solution statistics.
            // Use the given solver controls
            residualsIO solve
			(
				const dictionary& solverControls,
				PtrList<volScalarField>& sW,
				PtrList<volVectorField>& vW
			);

            //- Solve returning the solution statistics.
            // Solver controls read from fvSolution
            residualsIO solve
			(
				PtrList<volScalarField>& sW,
				PtrList<volVectorField>& vW
			);

            //- Solve returning the solution statistics.
            // Calculates only the increments without updating the variables
            // Used only for delta-form
            // Solver controls read from fvSolution
            residualsIO solveForIncr
			(
				PtrList<volScalarField>& sW,
				PtrList<volVectorField>& vW,
				PtrList<volScalarField>& dsW,
				PtrList<volVectorField>& dvW
			);

        friend Ostream& operator<<(Ostream& os, const coupledMatrix& cm)
        {
            os << cm.dSByS_ << endl;
            os << cm.dSByV_ << endl;
            os << cm.dVByS_ << endl;
            os << cm.dVByV_ << endl;

            return os;
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
