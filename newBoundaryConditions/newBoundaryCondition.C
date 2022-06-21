/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDist.H"
#include "wordRe.H"


int main(int argc, char *argv[])
{
    // timeSelector::addOptions();
    
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
	// runTime and mesh are instances of objects (or classes).

    volVectorField U_orig // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            mesh.time().path()/"0.orig",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField wallDistance_ = wallDist(mesh).y();


    // forAll(mesh.boundaryMesh(), patchI)
    // {
    //     Info << "patch size (faces): " << mesh.boundary()[patchI].patch().size() << endl;

    //     forAll(mesh.boundary()[patchI].patch(), patchFaceI)
    //     {
    //         label adjcellIdD = mesh.boundary()[patchI].patch().faceCells()[patchFaceI];
    //         Info << "Patch " << patchI << " has its face " << patchFaceI << " adjacent to cell "
    //              << adjcellIdD << ". The distance of the cell to wall is: "
    //              << wallDistance_[adjcellIdD] << endl;
    //     }        
    // }

    // Get index of patch
    wordRe matchPatch("inlet.*");
    matchPatch.compile ();

    List<label> matchPatchID(mesh.boundaryMesh().size());
    label countNumber(0);

    // wordList testL ({"inlet_1", "inlet_2", "inllet_3", "inle"});
    // forAll(testL, ID)
    // {
    //     Info << matchPatch.match(testL[ID]) << endl;
    // }

    forAll(mesh.boundaryMesh(), patchI)
    {
        // Info << matchPatch.match(mesh.boundary()[patchI].name()) << endl;
        // Info << "patch name is: " << mesh.boundary()[patchI].name() << endl;
        if(matchPatch.match(mesh.boundary()[patchI].name()))
        {
            // Info << "match patch name is: " << mesh.boundary()[patchI].name() << endl;
            // matchPatchID.append(patchI);
            matchPatchID[countNumber] = patchI;
            countNumber += 1;
        }
    }

    matchPatchID.resize(countNumber);

    // Info << matchPatchID << endl;

    // label patchI_ = mesh.boundaryMesh().findPatchID("inlet");

    // label patchI_ (0);

    // Info << patchI_ << endl;

    // fvPatchField<vector>& inletUorig = U_orig.boundaryFieldRef()[patchI_];

    // // loop over all hub faces
    // forAll(inletUorig, faceI)
    // {
    //     label adjcellID = mesh.boundary()[patchI_].patch().faceCells()[faceI];
    //     vector uValue(0, 0, 0.5*(pow(wallDistance_[adjcellID]/max(wallDistance_).value(), 4)));
    //     inletUorig[faceI] = uValue;
    // }

    Info << "max(wallDistance_).value(): " << max(wallDistance_).value() << endl;

    forAll(matchPatchID, matchI)
    {
        label patchI_ (matchPatchID[matchI]);
        // label patchI_ (0);

        Info << patchI_ << endl;

        fvPatchField<vector>& inletUorig = U_orig.boundaryFieldRef()[patchI_];

        // loop over all hub faces
        forAll(inletUorig, faceI)
        {
            label adjcellID = mesh.boundary()[patchI_].patch().faceCells()[faceI];
            // if (wallDistance_[adjcellID] > 5.0e-4)
            // {
            //     vector uValue(0, 0, 0.841*(pow(wallDistance_[adjcellID]/max(wallDistance_).value(), 2)));
            //     inletUorig[faceI] = uValue;
            // }            

            if (0.841*(pow(wallDistance_[adjcellID]/max(wallDistance_).value(), 2)) > 0.1)
            {
                vector uValue(0, 0, 0.841*(pow(wallDistance_[adjcellID]/max(wallDistance_).value(), 2)));
                inletUorig[faceI] = uValue;
            }
            else
            {
                vector uValue(0, 0, 0.1);
                inletUorig[faceI] = uValue;
            }
        }
    }
    
    U_orig.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
