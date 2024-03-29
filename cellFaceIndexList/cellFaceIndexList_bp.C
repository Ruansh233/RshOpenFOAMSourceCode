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
        "Specify the folder to read the center list"
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"   

    fileName dataPath (args.getOrDefault<word>("fileDir", "centerList"));
    
    // find the global cell ID of each cell in each cellZone
    List<vector> cellCenterList;

    autoPtr<IFstream> inputFilePtr;
    inputFilePtr.reset(new IFstream(runTime.globalPath()/dataPath/"cellCenterList"));
    inputFilePtr() >> cellCenterList;

    if (mesh.cellZones().size() == 0)
    {
        FatalErrorIn("cellZoneFieldValue")
            << "There is no cellZone in this mesh"
            << exit(FatalError);

    }

    List<List<label>> globalcellIndexList(mesh.cellZones().size());
    vector zoneCenter (0, 0, 0);

    
    // // test code for 1 zone
    // label zoneI = 0;
    // forAll(mesh.cellZones()[zoneI], cellI)
    // {
    //     label cellID = mesh.cellZones()[zoneI][cellI];
    //     zoneCenter += mesh.C()[cellID];
    // }
    // zoneCenter = zoneCenter/mesh.cellZones()[zoneI].size();

    // Info<< "zone_" << zoneI << " center: " << zoneCenter << endl;

    // globalcellIndexList[zoneI].resize(mesh.cellZones()[zoneI].size());
    // List<label> existCellID;

    // forAll(mesh.cellZones()[zoneI], cellI)
    // {     
    //     label cellID = mesh.cellZones()[zoneI][cellI];
        
    //     vector distanceToCenter (mesh.C()[cellID] - zoneCenter);
        
    //     forAll(cellCenterList, localI)
    //     {
    //         if (mag(distanceToCenter - cellCenterList[localI]) < 1.0e-8)
    //         {
    //             Info<< localI << " found" << endl;
    //             globalcellIndexList[zoneI][localI] = cellID;
    //             existCellID.append(localI);
    //             break;
    //         }        
    //     }
    // }


    forAll(mesh.cellZones(), zoneI)
    {
        vector zoneCenter (0, 0, 0);
        
        forAll(mesh.cellZones()[zoneI], cellI)
        {
            label cellID = mesh.cellZones()[zoneI][cellI];
            zoneCenter += mesh.C()[cellID];
        }
        zoneCenter = zoneCenter/mesh.cellZones()[zoneI].size();

        Info<< "zone_" << zoneI << " center: " << zoneCenter << endl;

        globalcellIndexList[zoneI].resize(mesh.cellZones()[zoneI].size());

        forAll(mesh.cellZones()[zoneI], cellI)
        {     
            label cellID = mesh.cellZones()[zoneI][cellI];
            
            vector distanceToCenter (mesh.C()[cellID] - zoneCenter);
            
            forAll(cellCenterList, localI)
            {
                if (mag(distanceToCenter - cellCenterList[localI]) < 1.0e-8)
                {
                    globalcellIndexList[zoneI][localI] = cellID;
                    break;
                }        
            }
        }
    }

    // output the cellZone List
    IOList<List<label>> globalcellIndexListIO
    (
        IOobject
        (
            "globalcellIndexList",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        globalcellIndexList
    );
    globalcellIndexListIO.write();  


    // find the global cell ID of each cell in each cellZone
    List <List<vector>> faceCenterList(mesh.boundary().size());

    autoPtr<IFstream> inputFilePtr;
    inputFilePtr.reset(new IFstream(runTime.globalPath()/dataPath/"cellCenterList"));
    inputFilePtr() >> cellCenterList;

    if (mesh.cellZones().size() == 0)
    {
        FatalErrorIn("cellZoneFieldValue")
            << "There is no cellZone in this mesh"
            << exit(FatalError);

    }

    List<List<label>> globalcellIndexList(mesh.cellZones().size());
    vector zoneCenter (0, 0, 0);




}


// ************************************************************************* //
