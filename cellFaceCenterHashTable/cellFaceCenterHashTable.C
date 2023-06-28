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

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "fileDir",
        "word",
        "Specify the folder to store the hash table"
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"   

    fileName dataPath (args.getOrDefault<word>("fileDir", "HashTable"));
    // create the SVD folder if it is not exist.
    if(!isDir(runTime.globalPath()/dataPath))
        mkDir(runTime.globalPath()/dataPath);
    
    const vector subdomainCenter (gAverage(mesh.C()));

    // the hash table of cell centers, label is the key（cell ID）, vector is the value(cell center coordinate)
    HashTable<label, vector> cellCenterTable;
    
    forAll(mesh.cells(), cellI)
    {
        vector distanceToCenter (mesh.C()[cellI] - subdomainCenter);

        if(mag(distanceToCenter.x()) < 1e-8)
            distanceToCenter.x() = 0;
        if(mag(distanceToCenter.y()) < 1e-8)
            distanceToCenter.y() = 0;
        if(mag(distanceToCenter.z()) < 1e-8)
            distanceToCenter.z() = 0;

        cellCenterTable.insert(distanceToCenter, cellI);
    }

    List< HashTable<label, vector> > gatheredcellCenterTable(Pstream::nProcs());
    gatheredcellCenterTable[Pstream::myProcNo()] = cellCenterTable;
    Pstream::gatherList(gatheredcellCenterTable);
    Pstream::scatterList(gatheredcellCenterTable);


    if (Pstream::master())
    {
        HashTable<label, vector> globalcellCenterTable;

        forAll(gatheredcellCenterTable, processorI)
        {
            globalcellCenterTable += gatheredcellCenterTable[processorI];
        }
        
        autoPtr<OFstream> outputFilePtr;
        outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"cellHashTable"));

        outputFilePtr() << globalcellCenterTable;
    }


    // surface hashtable
    forAll(mesh.boundary(), patchI)
    {
        word patchType = mesh.boundary()[patchI].type();
        
        wordRe patchMatch("processor.*");
        patchMatch.compile();

        if (!patchMatch.match(patchType))
        {
            Info<< "patchType_" << patchI << ": " << patchType << endl;
            
            HashTable<label, vector> faceCenterTable;

            forAll(mesh.boundary()[patchI].Cf(), faceI)
            {
                const vector distanceToCenter (mesh.boundary()[patchI].Cf()[faceI] - subdomainCenter);
                vector truncFaceCenter (trunc(distanceToCenter.x()*1e8),
                                        trunc(distanceToCenter.y()*1e8),
                                        trunc(distanceToCenter.z()*1e8));
                
                label faceID = mesh.boundary()[patchI].start() + faceI;
                faceCenterTable.insert(truncFaceCenter, faceID);
            }

            List< HashTable<label, vector> > gatheredfaceCenterTable(Pstream::nProcs());
            gatheredfaceCenterTable[Pstream::myProcNo()] = faceCenterTable;
            Pstream::gatherList(gatheredfaceCenterTable);
            Pstream::scatterList(gatheredfaceCenterTable);


            if (Pstream::master())
            {
                HashTable<label, vector> globalfaceCenterTable;

                forAll(gatheredfaceCenterTable, processorI)
                {
                    globalfaceCenterTable += gatheredfaceCenterTable[processorI];
                }
                
                autoPtr<OFstream> outputFilePtr;
                outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"faceHashTable_"+name(patchI)));

                outputFilePtr() << globalfaceCenterTable;
            }
        }
    }
}


// ************************************************************************* //
