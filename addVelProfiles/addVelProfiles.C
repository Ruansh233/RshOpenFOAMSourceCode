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

    // read the uniform ratio
    scalar a1(args.get<scalar>("a1"));
    scalar a2(args.get<scalar>("a2"));

    U = a1*U_1 + a2*U_2;
    
    U.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
