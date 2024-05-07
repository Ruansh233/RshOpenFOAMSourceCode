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
#include "IFstream.H"

int main(int argc, char *argv[])
{   
    argList::addOption
    (
        "patchNames",
        "word",
        "Specify the patchName for output, e.g. '(\"inlet\" \"outlet\")'"
    );

    argList::addBoolOption
    (
        "allPoints",
        "write all mesh points to file"
    );

    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    List<fileName> patchNames;    
    const bool patchBool (args.readListIfPresent <fileName> ("patchNames", patchNames));

    // output the pointField of patches
    if(patchBool)
    {
        forAll(mesh.boundaryMesh(), patchI)
        {
            if(patchNames.found(mesh.boundaryMesh()[patchI].name()))
            {
                Info<< "patch Name = " << mesh.boundaryMesh()[patchI].name() << endl;
                IOList<vector> patchPointIO
                (
                    IOobject
                    (
                        mesh.boundaryMesh()[patchI].name(),
                        runTime.caseConstant(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh.boundaryMesh()[patchI].localPoints()
                );
                patchPointIO.write();
            }
        }
    }
    else
    {
        Info<< "PatchNames is not specified, output all patches" << endl;
        forAll(mesh.boundaryMesh(), patchI)
        {
            Info<< "patch Name = " << mesh.boundaryMesh()[patchI].name() << endl;
            IOList<vector> patchPointIO
            (
                IOobject
                (
                    mesh.boundaryMesh()[patchI].name(),
                    runTime.caseConstant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh.boundaryMesh()[patchI].localPoints()
            );
            patchPointIO.write();
        }
    }

    if (args.found("allPoints"))
    {
        Info<< "Output all points" << endl;
        IOList<vector> allPointIO
        (
            IOobject
            (
                "points",
                runTime.caseConstant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh.points()
        );
        allPointIO.write();
    }
    

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
