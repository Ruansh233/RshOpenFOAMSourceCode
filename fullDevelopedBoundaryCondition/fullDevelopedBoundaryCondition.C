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

    volVectorField U_fd // note that velocity is a vector field
    (
        IOobject
        (
            "U_fd",
            mesh.time().path()/"0.orig",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Info << U_fd << endl;

    // Info << "mesh.time().path(): " << mesh.time().path()/"0" << endl;


    // Info << "runTime.timeName(): " << runTime.timeName() << endl
    //      << "mesh.time().timeName(): " << mesh.time().timeName() << endl
    //      << "startTime: " << runTime.startTime() << endl
    //      << "startTime: " << mesh.time().startTime() << endl
    //      << "path: " << mesh.time().path() << endl;

    // the output of forward Info
    // runTime.timeName(): 5, this is the startTime of the case
    // mesh.time().timeName(): 5, this is equal to runTime.timeName()
    // startTime: startTime [0 0 1 0 0 0 0] 5
    // startTime: startTime [0 0 1 0 0 0 0] 5
    // path: "/home/shenhui_ruan/gitFolder/RshOpenFOAMSourceCode/fullDevelopedBoundaryCondition/testCase"
        
    // volScalarField wallDistance_ = wallDist(mesh).y();

    // Get index of patch

    wordRe matchInlet("Block.*_bottom");
    wordRe matchOutlet("Block.*_top");


    // wordRe matchInlet("inlet");    
    matchInlet.compile ();

    // wordRe matchOutlet("outlet");    
    matchOutlet.compile ();

    List<label> matchInletID(mesh.boundaryMesh().size());
    label countNumberI(0);
    List<label> matchOutletID(mesh.boundaryMesh().size());
    label countNumberO(0);

    // wordList testL ({"inlet_1", "inlet_2", "inllet_3", "inle"});
    // forAll(testL, ID)
    // {
    //     Info << matchInlet.match(testL[ID]) << endl;
    // }

    forAll(mesh.boundaryMesh(), patchI)
    {
        // Info << matchInlet.match(mesh.boundary()[patchI].name()) << endl;
        // Info << "patch name is: " << mesh.boundary()[patchI].name() << endl;
        if(matchInlet.match(mesh.boundary()[patchI].name()))
        {
            // Info << "match patch name is: " << mesh.boundary()[patchI].name() << endl;
            // matchInletID.append(patchI);
            matchInletID[countNumberI] = patchI;
            countNumberI += 1;
        }

        if(matchOutlet.match(mesh.boundary()[patchI].name()))
        {
            // Info << "match patch name is: " << mesh.boundary()[patchI].name() << endl;
            // matchInletID.append(patchI);
            matchOutletID[countNumberO] = patchI;
            countNumberO += 1;
        }
    }

    matchInletID.resize(countNumberI);
    matchOutletID.resize(countNumberO);


    // Info << matchInletID << endl;

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

    forAll(matchInletID, matchI)
    {
        label patchI_ (matchInletID[matchI]);
        label patchO_ (matchOutletID[matchI]);
        // label patchI_ (0);

        Info << "inlet patch No: " << patchI_ << ". it's name is " << mesh.boundary()[patchI_].name() << endl
             << "outlet patch No: " << patchO_ << ". it's name is " << mesh.boundary()[patchO_].name() << endl;

        fvPatchField<vector>& inletU = U.boundaryFieldRef()[patchI_];
        fvPatchField<vector> inletU_fd = U_fd.boundaryField()[patchO_];

        // loop over all hub faces
        forAll(inletU, faceI)
        {
            inletU[faceI] = inletU_fd[faceI];
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
    
    U.write();

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
