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

// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone

#include "fvCFD.H"
#include "wordRe.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // patch and owner cellzone list
    List<List<word>> zonePatchName(mesh.cellZones().size());
    List<List<label>> zonePatchID(mesh.cellZones().size());

    Info<< "Start to find zones and related patches ID." << endl;
    
    // The stupid loop to find the zone and connected patch list
    forAll(mesh.cellZones(), zoneI)
    {
        word cellzoneName (mesh.cellZones()[zoneI].name());
        cellzoneName.removeEnd("Zone");
        wordRe matchZoneName("[a-zA-Z]*_?"+cellzoneName+".*");
        // wordRe matchZoneName(".*"+cellzoneName+".*");
        matchZoneName.compile();

        forAll(mesh.boundaryMesh(), patchI)
        {
            if(matchZoneName.match(mesh.boundary()[patchI].name()))
            {
                zonePatchID[zoneI].append(patchI);
                zonePatchName[zoneI].append(mesh.boundary()[patchI].name());
            }
        }
    }
    Info<< "Successing in finding zones and related patches ID." << endl;


    // IOList for list write
    IOList<List<label>> cellZonePatchIO
    (
        IOobject
        (
            "cellZonePatch",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        zonePatchID
    );

    cellZonePatchIO.write();    

    // IOList for list write
    IOList<List<word>> cellZonePatchNameIO
    (
        IOobject
        (
            "cellZonePatchName",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        zonePatchName
    );

    cellZonePatchNameIO.write();  

    // // test list read steps
    // // IOList for list input and output
    // IOList<List<label>> cellZonePatchIOTest
    // (
    //     IOobject
    //     (
    //         "cellZonePatch",
    //         runTime.time().system(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     )
    // );
    // Info<< cellZonePatchIOTest << endl; 

    Info << "\nEnd\n" << endl;

    return 0;

}


// ************************************************************************* //
