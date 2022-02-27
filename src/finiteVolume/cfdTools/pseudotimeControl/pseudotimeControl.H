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

#ifndef pseudotimeControl_H
#define pseudotimeControl_H

#include "solutionControl.H"
#include "residualsIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pseudotimeControl Declaration
\*---------------------------------------------------------------------------*/

class pseudotimeControl
:
    public solutionControl

{
protected:

    // Protected data

        // Solution controls

            //- Steady state flag
            bool steadyState_;

            //- Maximum number of outer (pseudotime) correctors
            label nCorrOuter_;

            //- Minimum number of outer (pseudotime) correctors
            label nCorrOuterMin_;

            //- Flag to indicate whether to only solve turbulence on final iter
            bool turbOnFinalIterOnly_;

            //- Converged flag
            bool converged_;

            //- First iteration or first iteration after restart (no residual computed yet)
            bool firstIteration_;

        // Residuals

            scalar residualTols_;
            scalar residualTolsRel_;

            residualsIO residuals_;
            residualsIO initResiduals_;

            // For persistent storage of initial residuals
            IOdictionary state_;

    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual bool read();

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied();

        //- Disallow default bitwise copy construct
        pseudotimeControl(const pseudotimeControl&);

        //- Disallow default bitwise assignment
        void operator=(const pseudotimeControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("pseudotimeControl");


    // Constructors

        //- Construct from mesh
        pseudotimeControl
        (
            fvMesh& mesh,
            const bool steadyState,
            const label nScalars,
            const label nVectors
        );


    //- Destructor
    virtual ~pseudotimeControl();


    // Member Functions

        // Access

            //- Maximum number of pseudotime correctors
            inline label nCorrOuter() const;

            //- Initial residual
			inline residualsIO& residuals();


        // Solution control

            //- pseudotime loop
            virtual bool loop();

            //- Helper function to identify when to store the intial residuals
            inline bool storeInitialResiduals() const;

            //- Helper function to identify first outer (pseudotime) iteration
            inline bool firstIter() const;

            //- Helper function to identify final outer (pseudotime) iteration
            inline bool finalIter() const;

            //- Helper function to identify whether to solve for turbulence
            inline bool turbCorr() const;

            //- Helper function to store residual (initial solver residual)
            inline void setResidual(const residualsIO& res);

            void setCorr(const label corr)
            {
                corr_ = corr;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pseudotimeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
