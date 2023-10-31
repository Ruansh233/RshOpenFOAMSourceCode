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
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"


    // -----------------------------------------------------------
    // --------- read cvSubMatchDict -----------------------------
    // -----------------------------------------------------------

    const word dictName("cvSubMatchDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;
    IOdictionary cvSubMatchDict(dictIO);

    // -----------------------------------------------------------
    // ------ find indices for control volume --------------------
    // -----------------------------------------------------------

    wordRe cvNameRex (cvSubMatchDict.getOrDefault<word>("cvNameRex", "cv.*_Zone"));
    cvNameRex.compile();
    const List<label> cvCellZones(mesh.cellZones().indices(cvNameRex));
    Info<< "cvCellZones: " << cvCellZones << endl;



    // -----------------------------------------------------------
    // --------- create new object for reference cases -----------
    // -----------------------------------------------------------
    // time object for reference cases
    Foam::Time runTimeTest
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        // refCaseName is ambigous if it type is word, so the fileName used
        refCaseName,
        "system",
        "constant"
    );

    // create new mesh object for reference cases
    fvMesh  refElementMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            args.rootPath()/refCaseName/"constant",
            runTimeTest,
            IOobject::MUST_READ
        ),
        false
    );

    wordRe cellZoneNameRex ("block.*_Zone");
    cellZoneNameRex.compile();

    const List<label> selectedCellZones(mesh.cellZones().indices(cellZoneNameRex));
    Info<< "selectedCellZones: " << selectedCellZones << endl;

    // -----------------------------------------------------------
    // --------- assign cell value for scalar field --------------
    // -----------------------------------------------------------

    forAll(scalarFields, i)
    {
        volScalarField scalarFieldRead
        (
            IOobject
            (
                scalarFields[i],
                runTimeTest.timeName(),
                refElementMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            refElementMesh
        );
        
        volScalarField scalarFieldWrite
        (
            IOobject
            (
                scalarFields[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
        );

        forAll(selectedCellZones, zoneI)
        {
            forAll(mesh.cellZones()[selectedCellZones[zoneI]], cellI)
            {
                label cell = mesh.cellZones()[selectedCellZones[zoneI]][cellI];
                scalarFieldWrite[cell] = scalarFieldRead[cellI];
            }
        }

        scalarFieldWrite.write();
    }

    forAll(vectorFields, i)
    {
        volVectorField vectorFieldRead
        (
            IOobject
            (
                vectorFields[i],
                runTimeTest.timeName(),
                refElementMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            refElementMesh
        );
        
        volVectorField vectorFieldWrite
        (
            IOobject
            (
                vectorFields[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
        );

        forAll(selectedCellZones, zoneI)
        {
            forAll(mesh.cellZones()[selectedCellZones[zoneI]], cellI)
            {
                label cell = mesh.cellZones()[selectedCellZones[zoneI]][cellI];
                vectorFieldWrite[cell].x() = vectorFieldRead[cellI].x();
                vectorFieldWrite[cell].y() = vectorFieldRead[cellI].y();
                vectorFieldWrite[cell].z() = vectorFieldRead[cellI].z();
            }
        }
        vectorFieldWrite.write();
    }

    forAll(tensorFields, i)
    {
        volTensorField tensorFieldRead
        (
            IOobject
            (
                tensorFields[i],
                runTimeTest.timeName(),
                refElementMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            refElementMesh
        );
        
        volTensorField tensorFieldWrite
        (
            IOobject
            (
                tensorFields[i],
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
        );

        forAll(selectedCellZones, zoneI)
        {
            forAll(mesh.cellZones()[selectedCellZones[zoneI]], cellI)
            {
                label cell = mesh.cellZones()[selectedCellZones[zoneI]][cellI];
                tensorFieldWrite[cell].xx() = tensorFieldRead[cellI].xx();
                tensorFieldWrite[cell].xy() = tensorFieldRead[cellI].xy();
                tensorFieldWrite[cell].xz() = tensorFieldRead[cellI].xz();
                tensorFieldWrite[cell].yx() = tensorFieldRead[cellI].yx();
                tensorFieldWrite[cell].yy() = tensorFieldRead[cellI].yy();
                tensorFieldWrite[cell].yz() = tensorFieldRead[cellI].yz();
                tensorFieldWrite[cell].zx() = tensorFieldRead[cellI].zx();
                tensorFieldWrite[cell].zy() = tensorFieldRead[cellI].zy();
                tensorFieldWrite[cell].zz() = tensorFieldRead[cellI].zz();
            }
        }
        tensorFieldWrite.write();
    }

}
    