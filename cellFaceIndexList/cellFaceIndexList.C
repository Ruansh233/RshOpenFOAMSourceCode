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

    argList::addBoolOption
    (
        "no-cell",
        "whether to store the list of cell centers or not"
    );

    argList::addBoolOption
    (
        "no-face",
        "whether to store the list of face centers or not"
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"   

    fileName dataPath (args.getOrDefault<word>("fileDir", "centerList"));
    
    List<vector> zoneCenterList(mesh.cellZones().size());
    vector zoneCenter (0, 0, 0);

    forAll(mesh.cellZones(), zoneI)
    {
        vector zoneCenter (0, 0, 0);
      
        forAll(mesh.cellZones()[zoneI], cellI)
        {
            label cellID = mesh.cellZones()[zoneI][cellI];
            zoneCenter += mesh.C()[cellID];
        }

        zoneCenterList[zoneI] = zoneCenter/mesh.cellZones()[zoneI].size();
    }


    if(!args.found("no-cell"))
    {
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

        forAll(mesh.cellZones(), zoneI)
        {
            Info<< "zone_" << zoneI << " center: " << zoneCenterList[zoneI] << endl;

            globalcellIndexList[zoneI].resize(mesh.cellZones()[zoneI].size());

            forAll(mesh.cellZones()[zoneI], cellI)
            {     
                label cellID = mesh.cellZones()[zoneI][cellI];
                
                vector distanceToCenter (mesh.C()[cellID] - zoneCenterList[zoneI]);
                
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
        List< List<List<label>> > gatheredcellIndexList(Pstream::nProcs());
        gatheredcellIndexList[Pstream::myProcNo()] = globalcellIndexList;
        Pstream::gatherList(gatheredcellIndexList);
        Pstream::scatterList(gatheredcellIndexList);

        if (Pstream::master())
        {
            List<List<label>> gatheredglobalcellIndexList  = 
            ListListOps::combine< List<List<label>> > 
            (gatheredcellIndexList, accessOp< List<List<label>> >());
            
            autoPtr<OFstream> outputFilePtr;
            outputFilePtr.reset(new OFstream(runTime.caseConstant()/"globalcellIndexList"));

            outputFilePtr() << gatheredglobalcellIndexList;
        }
    }
    

    if(!args.found("no-face"))
    {
        // list of face centers
        List <List<vector>> faceCenterList(mesh.boundary().size());

        autoPtr<IFstream> inputFilePtr;
        inputFilePtr.reset(new IFstream(runTime.globalPath()/dataPath/"faceCenterList"));
        inputFilePtr() >> faceCenterList;

        if (mesh.faceZones().size() == 0)
        {
            FatalErrorIn("faceZoneFieldValue")
                << "There is no faceZone in this mesh"
                << exit(FatalError);
        }
        

        forAll(faceCenterList, patchI)
        {
            Info<< "patch_" << patchI << endl;

            List<List<label>> globalfaceIndexList(mesh.faceZones().size());
            
            forAll(mesh.faceZones(), zoneI)
            {           
                globalfaceIndexList[zoneI].resize(faceCenterList[patchI].size());

                forAll(mesh.faceZones()[zoneI], faceI)
                {     
                    label faceID = mesh.faceZones()[zoneI][faceI];
                    
                    vector distanceToCenter (mesh.Cf()[faceID] - zoneCenterList[zoneI]);
                    
                    forAll(faceCenterList[patchI], localI)
                    {
                        if (mag(distanceToCenter - faceCenterList[patchI][localI]) < 1.0e-8)
                        {
                            globalfaceIndexList[zoneI][localI] = faceID;
                            break;
                        }        
                    }
                }
            }

            // output the faceZone List
            List< List<List<label>> > gatheredfaceIndexList(Pstream::nProcs());
            gatheredfaceIndexList[Pstream::myProcNo()] = globalfaceIndexList;
            Pstream::gatherList(gatheredfaceIndexList);
            Pstream::scatterList(gatheredfaceIndexList);


            if (Pstream::master())
            {
                List<List<label>> gatheredglobalfaceIndexList  = 
                ListListOps::combine< List<List<label>> > 
                (gatheredfaceIndexList, accessOp< List<List<label>> >());
                
                autoPtr<OFstream> outputFilePtr;
                outputFilePtr.reset(new OFstream(runTime.caseConstant()/"gatheredfaceIndexList"));

                outputFilePtr() << gatheredglobalfaceIndexList;
            }
        }
    }
}


// ************************************************************************* //
