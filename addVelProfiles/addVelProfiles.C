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
#include "IFstream.H"
#include "stringOps.H"


int main(int argc, char *argv[])
{
    argList::addOption //
    (
        "a1",
        "scalar",
        "the ratio for uniform distribution."
    );

    argList::addOption //
    (
        "a2",
        "scalar",
        "the ratio for non-uniform distribution."
    );

    argList::addOption //
    (
        "DirIndex",
        "labelList",
        "index of Dirichlet boundary patches, e.g., '(0)' or '(0 1)'."
    );

    
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
	// runTime and mesh are instances of objects (or classes).

    volVectorField U_1 // note that velocity is a vector field
    (
        IOobject
        (
            "U_1",
            mesh.time().path()/"0.orig",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField U_2 // note that velocity is a vector field
    (
        IOobject
        (
            "U_2",
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
            mesh.time().path()/"0",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    labelList DirIndex (args.getList<label>("DirIndex"));

    // read the uniform ratio
    scalar a1(args.get<scalar>("a1"));
    scalar a2(args.get<scalar>("a2"));

    // print the two ratios
    Info<< "The ratios are: "
        << "a1 = " << a1 << ", "
        << "a2 = " << a2 << endl;

    // Assign cell values
    U = a1*U_1 + a2*U_2;

    // Assign Dirichlet BCs
    forAll(DirIndex, indexI)
    {
        label patchI_ (DirIndex[indexI]);

        // print the Dirichlet BCs
        Info<< "Dirichlet BC No." << patchI_ << " is: " 
            << mesh.boundary()[patchI_].name() << endl;

        fvPatchField<vector>& inletUorig = U.boundaryFieldRef()[patchI_];
        const vectorField& U1_patch = U_1.boundaryField()[patchI_];
        const vectorField& U2_patch = U_2.boundaryField()[patchI_];

        // loop over all hub faces
        forAll(mesh.boundary()[patchI_].Cf(), faceI)
        {
            inletUorig[faceI] = a1*U1_patch[faceI] + a2*U2_patch[faceI];
        }
    }
    
    U.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
