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
        "Specify the folder to store the list of cell centers and face centers."
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"   

    fileName dataPath (args.getOrDefault<word>("fileDir", "centerList"));
    // create the SVD folder if it is not exist.
    if(!isDir(runTime.globalPath()/dataPath))
        mkDir(runTime.globalPath()/dataPath);
    
    const vector subdomainCenter (gAverage(mesh.C()));

    // the list of cell centers
    List<vector> cellCenterList;
    
    forAll(mesh.cells(), cellI)
    {
        vector distanceToCenter (mesh.C()[cellI] - subdomainCenter);

        cellCenterList.append(distanceToCenter);
    }

    List< List<vector> > gatheredcellCenterList(Pstream::nProcs());
    gatheredcellCenterList[Pstream::myProcNo()] = cellCenterList;
    Pstream::gatherList(gatheredcellCenterList);
    Pstream::scatterList(gatheredcellCenterList);


    if (Pstream::master())
    {
        List<vector> globalcellCenterList  = 
        ListListOps::combine<List<vector>> 
        (gatheredcellCenterList, accessOp<List<vector> >());
        
        autoPtr<OFstream> outputFilePtr;
        outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"cellCenterList"));

        outputFilePtr() << globalcellCenterList;
    }


    // list of face centers
    List <List<vector>> faceCenterList(mesh.boundary().size());
    
    forAll(mesh.boundary(), patchI)
    {
        faceCenterList[patchI].resize(mesh.boundary()[patchI].size());

        word patchType = mesh.boundary()[patchI].type();
        
        wordRe patchMatch("processor.*");
        patchMatch.compile();

        if (!patchMatch.match(patchType))
        {
            Info<< "patchType_" << patchI << ": " << patchType << endl;

            forAll(mesh.boundary()[patchI].Cf(), faceI)
            {
                const vector distanceToCenter (mesh.boundary()[patchI].Cf()[faceI] - subdomainCenter);
                
                faceCenterList[patchI][faceI] = distanceToCenter;
            }
        }
    }

    List< List <List<vector>> > gatheredfaceCenterList(Pstream::nProcs());
    gatheredfaceCenterList[Pstream::myProcNo()] = faceCenterList;
    Pstream::gatherList(gatheredfaceCenterList);
    Pstream::scatterList(gatheredfaceCenterList);

    if (Pstream::master())
    {
        List <List<vector>> globalfaceCenterList  = 
        ListListOps::combine<List <List<vector>> > 
        (gatheredfaceCenterList, accessOp<List <List<vector>> >());
        
        autoPtr<OFstream> outputFilePtr;
        outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"faceCenterList"));

        outputFilePtr() << globalfaceCenterList;
    }
}


// ************************************************************************* //
