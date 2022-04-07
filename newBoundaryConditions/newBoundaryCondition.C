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
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

	// runTime and mesh are instances of objects (or classes).
	// If you are not familiar with what a class or object is, it is HIGHLY RECOMMENDED you visit this
	// website and only come back once you've read everything about classes, inheritance and polymorphism:
	// http://www.cplusplus.com/doc/tutorial/classes/
	// Note how the next lines call functions .timeName(), .C() and .Cf() implemented in the objects.
	// It is also important to realise that mesh.C() and .Cf() return vector fields denoting centres of each
    // cell and internal face.
	// Calling the mesh.C().size() method therefore yields the total size of the mesh.

    
	// Info << "Hello there, the most recent time folder found is " << runTime.timeName() << nl
	// 	 << "The mesh has " << mesh.C().size() << " cells and " << mesh.Cf().size()
    //      << " internal faces in it. " << nl << endl;

    // // It's possible to iterate over every cell in a standard C++ for loop
    // for (label cellI = 0; cellI < mesh.C().size(); cellI++)
    //     if (cellI%20 == 0) // only show every twentieth cell not to spam the screen too much
    //         Info << "Cell " << cellI << " with centre at " << mesh.C()[cellI] << endl;
    // Info << endl; // spacer

    // // Each cell is constructed of faces - these may either be internal or constitute a
    // // boundary, or a patch in OpenFOAM terms; internal faces have an owner cell
    // // and a neighbour.
    // for (label faceI = 0; faceI < mesh.owner().size(); faceI++)
    //     if (faceI%40 == 0)
    //         Info << "Internal face " << faceI << " with centre at " << mesh.Cf()[faceI]
    //              << " with owner cell " << mesh.owner()[faceI]
    //              << " and neighbour " << mesh.neighbour()[faceI] << endl;
    // Info << endl;

    // // Boundary conditions may be accessed through the boundaryMesh object.
    // // In reality, each boundary face is also included in the constant/polyMesh/faces
    // // description. But, in that file, the internal faces are defined first.
    // // In addition, the constant/polyMesh/boundary file defines the starting faceI
    // // indices from which boundary face definitions start.
    // // OpenFOAM also provides a macro definition for for loops over all entries
    // // in a field or a list, which saves up on the amount of typing.
    // forAll(mesh.boundaryMesh(), patchI)
    //     Info << "Patch " << patchI << ": " << mesh.boundary()[patchI].name() << " with "
    //          << mesh.boundary()[patchI].Cf().size() << " faces. Starts at total face "
    //          << mesh.boundary()[patchI].start() << endl;
    // Info << endl;

    // // Faces adjacent to boundaries may be accessed as follows.
    // // Also, a useful thing to know about a face is its normal vector and face area.
    // label patchFaceI(0);
    // forAll(mesh.boundaryMesh(), patchI)
    //     Info << "Patch " << patchI << " has its face " << patchFaceI << " adjacent to cell "
    //          << mesh.boundary()[patchI].patch().faceCells()[patchFaceI]
    //          << ". It has normal vector " << mesh.boundary()[patchI].Sf()[patchFaceI]
    //          << " and surface area " << mag(mesh.boundary()[patchI].Sf()[patchFaceI])
    //          << endl;
    // Info << endl;

    // fileName veloctityIni (runTime.path()/"0.orig");

    // volVectorField U_orig // note that velocity is a vector field
    // (
    //     IOobject
    //     (
    //         "U",
    //         veloctityIni,
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh
    // );

    volVectorField U_orig // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTime.timeName(),
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
            vector uValue(0, 0, 1.0*(pow(wallDistance_[adjcellID]/max(wallDistance_).value(), 2)));
            inletUorig[faceI] = uValue;
        }
    }

    // label inletPatchID = mesh.boundaryMesh().findPatchID("inlet.*");

    // Get reference to boundary value
    // const fvPatchVectorField& faceCentreshub = mesh.Cf().boundaryField()[inletPatchID];
    // fvPatchField<vector>& inletUorig = U_orig.boundaryFieldRef()[inletPatchID];

    // Info << "patch " << inletPatchID << ": " << mesh.boundary()[inletPatchID].name() 
    //          <<" has " << mesh.boundary()[inletPatchID].patch().size() << " faces." << endl;

    // // loop over all hub faces
    // forAll(inletUorig, faceI)
    // {
    //     // // get coordinate for face centre
    //     // const vector& c = faceCentreshub[faceI];
    //     // vector p(0.5*(1+Foam::sin(40*M_PI*c[0]-M_PI/2)), 0, 0);

    //     // if (c[0]>0.025 &c[0]<0.075)
    //     //     p = vector(1, 0, 0);

    //     // movingWallU[faceI] = p;

    //     label adjcellID = mesh.boundary()[inletPatchID].patch().faceCells()[faceI];
    //     // vector uValue(0, 0, 0.5*(wallDistance_[adjcellID]/max(wallDistance_).value()));
    //     vector uValue(0, 0, 0.5*(pow(wallDistance_[adjcellID]/max(wallDistance_).value(), 4)));


    //     inletUorig[adjcellID] = uValue;
    // }

    // volVectorField U // note that velocity is a vector field
    // (
    //     IOobject
    //     (
    //         "U",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh,
    //     dimVelocity,
    //     U_orig
    // );
    
    U_orig.write();

    // // For internal faces, method .Sf() can be called directly on the mesh instance.
    // // Moreover, there is a shorthand method .magSf() which returns the surface area
    // // as a scalar.
    // // For internal faces, the normal vector points from the owner to the neighbour
    // // and the owner has a smaller cellI index than the neighbour. For boundary faces,
    // // the normals always point outside of the domain (they have "imaginary" neighbours
    // // which do not exist).

    // // It is possible to look at the points making up each face in more detail.
    // // First, we define a few shorthands by getting references to the respective
    // // objects in the mesh. These are defined as constants since we do not aim to
    // // alter the mesh in any way.
    // // NOTE: these lists refer to the physical definition of the mesh and thus
    // // include boundary faces. Use can be made of the mesh.boundary()[patchI].Cf().size()
    // // and mesh.boundary()[patchI].start() methods to check whether the face is internal
    // // or lies on a boundary.

    // const faceList& fcs = mesh.faces();
    // const List<point>& pts = mesh.points();
    // const List<point>& cents = mesh.faceCentres();

    // forAll(fcs,faceI)
    //     if (faceI%80==0)
    //     {
    //         if (faceI<mesh.Cf().size())
    //             Info << "Internal face ";
    //         else
    //         {
    //             forAll(mesh.boundary(),patchI)
    //                 if ((mesh.boundary()[patchI].start()<= faceI) &&
    //                     (faceI < mesh.boundary()[patchI].start()+mesh.boundary()[patchI].Cf().size()))
    //                 {
    //                     Info << "Face on patch " << patchI << ", faceI ";
    //                     break; // exit the forAll loop prematurely
    //                 }
    //         }

    //         Info << faceI << " with centre at " << cents[faceI]
    //              << " has " << fcs[faceI].size() << " vertices:";
    //         forAll(fcs[faceI],vertexI)
    //             // Note how fcs[faceI] holds the indices of points whose coordinates
    //             // are stored in the pts list.
    //             Info << " " << pts[fcs[faceI][vertexI]];
    //         Info << endl;
    //     }
    // Info << endl;

    // // In the original cavity tutorial, on which the test case is based,
    // // the frontAndBack boundary is defined as and "empty" type. This is a special
    // // BC case which may cause unexpected behaviour as its .Cf() field has size of 0.
    // // Type of a patch may be checked to avoid running into this problem if there
    // // is a substantial risk that an empty patch type will appear

    // label patchID(0);
    // const polyPatch& pp = mesh.boundaryMesh()[patchID];
    // if (isA<emptyPolyPatch>(pp))
    // {
    //     // patch patchID is of type "empty".
    //     Info << "You will not see this." << endl;
    // }

    // // Patches may also be retrieved from the mesh using their name. This could be
    // // useful if the user were to refer to a particular patch from a dictionary
    // // (like when you do when calculating forces on a particular patch).
    // word patchName("movingWall");
    // patchID = mesh.boundaryMesh().findPatchID(patchName);
    // Info << "Retrieved patch " << patchName << " at index " << patchID << " using its name only." << nl << endl;

    // Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
