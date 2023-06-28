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
#include "HashTable.H"
#include "IFstream.H"
#include "IOmanip.H"

int main(int argc, char *argv[])
{    
    argList::addOption
    (
        "fileDir",
        "word",
        "Specify the folder to read the hash table"
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"   

    fileName dataPath (args.getOrDefault<word>("fileDir", "HashTable"));
    HashTable<label, vector> cellHashTable;

    autoPtr<IFstream> inputFilePtr;
    inputFilePtr.reset(new IFstream(runTime.globalPath()/dataPath/"cellHashTable"));
    inputFilePtr() >> cellHashTable;

    if (mesh.cellZones().size() == 0)
    {
        FatalErrorIn("cellZoneFieldValue")
            << "There is no cellZone in this mesh"
            << exit(FatalError);

    }

    // List<List<label>> globalIndexList(mesh.cellZones().size());
    List<List<label>> globalIndexList;

    vector zoneCenter (0, 0, 0);
    label zoneI = 0;

    // Info<< "name: " << mesh.cellZones()[zoneI].name() << endl
    //     << "size: " << mesh.cellZones()[zoneI].size() << endl;

    forAll(mesh.cellZones()[zoneI], cellI)
    {
        label cellID = mesh.cellZones()[zoneI][cellI];
        zoneCenter += mesh.C()[cellID];
    }
    zoneCenter = zoneCenter/mesh.cellZones()[zoneI].size();

    Info<< "zone_" << zoneI << " center: " << zoneCenter << endl;

    List<label> globalCellIndex(mesh.cellZones()[zoneI].size());

    forAll(mesh.cellZones()[zoneI], cellI)
    {     
        label cellID = mesh.cellZones()[zoneI][cellI];
        
        vector distanceToCenter (mesh.C()[cellID] - zoneCenter);
        if(mag(distanceToCenter.x()) < 1e-8)
            distanceToCenter.x() = 0;
        if(mag(distanceToCenter.y()) < 1e-8)
            distanceToCenter.y() = 0;
        if(mag(distanceToCenter.z()) < 1e-8)
            distanceToCenter.z() = 0;
        
        Info<< "cell " << cellID << "cellcenter: " << mesh.C()[cellID]  << ", distanceToCenter: " << distanceToCenter << endl;     
        
        if (cellHashTable.found(distanceToCenter))
            globalCellIndex[cellI] = cellID;
        else
            globalCellIndex[cellI] = -1;
    }

    // globalIndexList[zoneI] = globalCellIndex;
    globalIndexList.append (globalCellIndex);


    // output the cellZone List
    IOList<List<label>> globalIndexListIO
    (
        IOobject
        (
            "globalIndexList",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        globalIndexList
    );
    globalIndexListIO.write();  
}


// ************************************************************************* //
