/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
    

Group
    

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "searchableSurface.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "listFolder",
        "name",
        "folder to save list of matching"
    );

    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative cvSubMatchDict"
    );

    argList::addOption
    (
        "scale toleration",
        "scalar",
        "toleration for scale the triSurfaceMesh, default is 1.0e-6"
    );
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"


    // -----------------------------------------------------------
    // --------- input & output initialization -------------------
    // -----------------------------------------------------------

    autoPtr<OFstream> outputFilePtr;
    fileName dataPath;
    scalar scaleToleration = args.getOrDefault<scalar>("scale toleration", 1.0e-6);
    
    // -----------------------------------------------------------
    // --------- read cvSubMatchDict -----------------------------
    // -----------------------------------------------------------

    const word dictName("cvSubMatchDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;
    IOdictionary cvSubMatchDict(dictIO);

    List<word> triSurfaceNames
    (
        cvSubMatchDict.lookup("triSurfaceNames")
    );

    
    // -----------------------------------------------------------
    // --------- list for triSurfaces ----------------------------
    // -----------------------------------------------------------
    PtrList<triSurfaceMesh> subTriSurfaces;
    forAll(triSurfaceNames, triSurfaceI)
    {
        subTriSurfaces.append
        (new triSurfaceMesh 
        (
            IOobject
            (
                triSurfaceNames[triSurfaceI],
                mesh.time().constant(),     
                "triSurface",              
                mesh.objectRegistry::db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ));
    }

    // PtrList<triSurfaceMesh> subTriSurfaces(triSurfaceNames.size());
    // forAll(triSurfaceNames, triSurfaceI)
    // {
    //     subTriSurfaces.set
    //     (
    //         triSurfaceI,
    //         new triSurfaceMesh
    //         (
    //             IOobject
    //             (
    //                 triSurfaceNames[triSurfaceI],
    //                 mesh.time().constant(),     
    //                 "triSurface",              
    //                 mesh.objectRegistry::db(),
    //                 IOobject::MUST_READ,
    //                 IOobject::NO_WRITE
    //             )
    //         )
    //     );
    // }

    

    // -----------------------------------------------------------
    // ------ find indices for control volume cell ---------------
    // -----------------------------------------------------------

    wordRe cvNameRex (cvSubMatchDict.getOrDefault<word>("cvNameRex", "cv.*_Zone"));
    cvNameRex.compile();
    const List<label> cvCellZones(mesh.cellZones().indices(cvNameRex));
    Info<< "cvCellZones: " << cvCellZones << endl;
    List<List<label>> cvCellZoneMatchID(cvCellZones.size());

    forAll(cvCellZones, cvI)
    {
        const label cvCellZone = cvCellZones[cvI];
        const labelList& cvCells = mesh.cellZones()[cvCellZone];
        const label cvCellZoneSize = cvCells.size();
        cvCellZoneMatchID[cvI].setSize(cvCellZoneSize);
        cvCellZoneMatchID[cvI] = -1;

        pointField cvCellCentres(cvCellZoneSize);
        forAll(cvCells, cellI)
        {
            cvCellCentres[cellI] = mesh.C()[cvCells[cellI]];
        }

        forAll(subTriSurfaces, surafaceI)
        {
            List<volumeType> volTypes;
            subTriSurfaces[surafaceI].scalePoints(1+scaleToleration);
            subTriSurfaces[surafaceI].getVolumeType(cvCellCentres, volTypes);

            forAll(volTypes, elemi)
            {
                if (volTypes[elemi] == volumeType::INSIDE)
                {
                    cvCellZoneMatchID[cvI][elemi] = surafaceI;
                }
            }
        }
    }

    dataPath = args.getOrDefault<word>("listFolder", "constant");
    outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"cvCellZoneMatchID"));

    outputFilePtr() << cvCellZoneMatchID;


    // -----------------------------------------------------------
    // --------- find indices for control volume face ------------
    // -----------------------------------------------------------

    wordRe cvFaceNameRex (cvSubMatchDict.getOrDefault<word>("cvFaceNameRex", "cv.*_face"));
    cvFaceNameRex.compile();
    const List<label> cvFaceZones(mesh.faceZones().indices(cvFaceNameRex));
    Info<< "cvFaceZones: " << cvFaceZones << endl;
    List<List<label>> cvFaceZoneMatchID(cvFaceZones.size());

    forAll(cvFaceZones, cvI)
    {
        const label cvFaceZone = cvFaceZones[cvI];
        const labelList& cvFaces = mesh.faceZones()[cvFaceZone];
        const label cvFaceZoneSize = cvFaces.size();
        cvFaceZoneMatchID[cvI].setSize(cvFaceZoneSize);
        cvFaceZoneMatchID[cvI] = -1;

        pointField cvFaceCentres(cvFaceZoneSize);
        forAll(cvFaces, faceI)
        {
            cvFaceCentres[faceI] = mesh.Cf()[cvFaces[faceI]];
        }

        forAll(subTriSurfaces, surafaceI)
        {
            List<volumeType> volTypes;
            subTriSurfaces[surafaceI].scalePoints(1+scaleToleration);
            subTriSurfaces[surafaceI].getVolumeType(cvFaceCentres, volTypes);

            forAll(volTypes, elemi)
            {
                if (volTypes[elemi] == volumeType::INSIDE)
                {
                    cvFaceZoneMatchID[cvI][elemi] = surafaceI;
                }
            }
        }
    }

    dataPath = args.getOrDefault<word>("listFolder", "constant");
    outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"cvFaceZoneMatchID"));

    outputFilePtr() << cvFaceZoneMatchID;


    Info<< "End" << endl;   
}
    