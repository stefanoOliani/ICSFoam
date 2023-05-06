/*---------------------------------------------------------------------------*\
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.

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
    Foam::MRFTranslatingZone

Description
    MRF zone definition based on cell zone and parameters
    obtained from a control dictionary constructed from the given stream.

SourceFiles
    MRFTranslatingZone.C

\*---------------------------------------------------------------------------*/

#ifndef MRFTranslatingZone_H
#define MRFTranslatingZone_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "mapPolyMesh.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class MRFTranslatingZone Declaration
\*---------------------------------------------------------------------------*/

class MRFTranslatingZone
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Name of the MRF region
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;

        //- MRF region active flag
        bool active_;

        //- Name of cell zone
        word cellZoneName_;

        //- Cell zone ID
        label cellZoneID_;

        const wordRes excludedPatchNames_;

        labelList excludedPatchLabels_;

        //- Internal faces that are part of MRF
        labelList internalFaces_;

        //- Outside faces (per patch) that move with the MRF
        labelListList includedFaces_;

        //- Excluded faces (per patch) that do not move with the MRF
        labelListList excludedFaces_;

        //- Velocity of the frame (rad/sec)
        autoPtr<Function1<vector>> frameVel_;


    // Private Member Functions

        //- Divide faces in frame according to patch
        void setMRFFaces();

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            FieldField<fvsPatchField, scalar>& phi
        ) const;

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            Field<scalar>& phi,
            const label patchi
        ) const;

        //- Make the given relative mass/vol flux absolute within the MRF region
        template<class RhoFieldType>
        void makeAbsoluteRhoFlux
        (
            const RhoFieldType& rho,
            surfaceScalarField& phi
        ) const;

        //- No copy construct
        MRFTranslatingZone(const MRFTranslatingZone&) = delete;

        //- No copy assignment
        void operator=(const MRFTranslatingZone&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("MRFTranslatingZone");


    // Constructors

        //- Construct from fvMesh
        MRFTranslatingZone
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName = word::null
        );

        //- Return clone
        autoPtr<MRFTranslatingZone> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Member Functions

        //- Return const access to the MRF region name
        inline const word& name() const;

        //- Return const access to the MRF active flag
        inline bool active() const;

        //- Return the current Omega vector
        vector frameVelocity() const;

        //- Update the mesh corresponding to given map
        void updateMesh(const mapPolyMesh& mpm)
        {
            // Only updates face addressing
            setMRFFaces();
        }

        //- Add omega within the MRF region
        void addFrameVelocity(volVectorField& UFrame) const;

        //- Make the given absolute velocity relative within the MRF region
        void makeRelative(volVectorField& U) const;

        //- Make the given absolute flux relative within the MRF region
        void makeRelative(surfaceScalarField& phi) const;

        //- Make the given absolute boundary flux relative
        //  within the MRF region
        void makeRelative(FieldField<fvsPatchField, scalar>& phi) const;

        //- Make the given absolute patch flux relative
        //  within the MRF region
        void makeRelative(Field<scalar>& phi, const label patchi) const;

        //- Make the given absolute mass-flux relative within the MRF region
        void makeRelative
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given relative velocity absolute within the MRF region
        void makeAbsolute(volVectorField& U) const;

        //- Make the given relative flux absolute within the MRF region
        void makeAbsolute(surfaceScalarField& phi) const;

        //- Make the given relative mass-flux absolute within the MRF region
        void makeAbsolute
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Compute the pseudo face velocity of the MRF region
	    void faceU(surfaceVectorField& zoneFaceU) const;

        //- Correct the boundary velocity for the rotation of the MRF region
        void correctBoundaryVelocity(volVectorField& U) const;

        //- Zero the MRF region of the given field
        template<class Type>
        void zero(GeometricField<Type, fvsPatchField, surfaceMesh>& phi) const;

        //- Update MRFTranslatingZone faces if the mesh topology changes
        void update();


    // I-O

        //- Write
        void writeData(Ostream& os) const;

        //- Read MRF dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MRFTranslatingZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MRFTranslatingZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
