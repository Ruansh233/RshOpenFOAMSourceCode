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

// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone

#include "fvCFD.H"
#include "wordRe.H"

int main(int argc, char *argv[])
{
    // Rsh, add dict select function
    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative cellZoneDict"
    );

    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // add dict Selector option
    const word dictName("cellZoneDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

    // create dictionary object and read the dictionary
    IOdictionary customDict(dictIO);
    #include "readDict.H"
    #include "createFields.H"


    word outputDir(customDict.getWord("outputDir"));
    // create the SVD folder if it is not exist.
    if(!isDir(runTime.globalPath()/outputDir))
        mkDir(runTime.globalPath()/outputDir);


    fileName snapshotFileName;
    autoPtr<OFstream> outputFilePtr;

    // p, pressure cell value
    forAll(typesName, typeI)
    {
        label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
        Info<< "Size of the snapshots, p" << name(typeI) << ", " 
            << zoneCells << ", " << cellZonesType[typeI].size() << endl;

        snapshotFileName = "p"+name(typeI);
        outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

        forAll(cellZonesType[typeI], zoneI)
        {                   
            forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
            {
                label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
                outputFilePtr() << p[cellN] << " ";
            }
            outputFilePtr() << endl;  
        } 
    }


    // vectorFields cell value
    forAll(typesName, typeI)
    {           
        label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
        Info<< "Size of the snapshots, U" << name(typeI) << ", "
            << zoneCells << ", " << cellZonesType[typeI].size() << endl;

        snapshotFileName = "U"+name(typeI);
        outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

        forAll(cellZonesType[typeI], zoneI)
        {
            forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
            {
                label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
                outputFilePtr() << U[cellN].x() << " "; 
            }

            forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
            {
                label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
                outputFilePtr() << U[cellN].y() << " "; 
            }

            forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
            {
                label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
                outputFilePtr() << U[cellN].z() << " "; 
            }
            
            outputFilePtr() << endl;   
        } 
    }

    
    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
