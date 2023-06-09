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
#include "emptyPolyPatch.H"
#include "IFstream.H"
#include "stringOps.H"

int main(int argc, char *argv[])
{
    // Rsh, add dict select function
    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative reconstructDict"
    );

    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    // add dict Selector option
    const word dictName("reconstructDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

    // create dictionary object and read the dictionary
    IOdictionary customDict(dictIO);
    // read data from dictionary
    #include "readDictList.H"

    // write scalar field field
    forAll(scalarZoneFieldName, nameNo)
    {
        Info << "scalarZoneFieldName: " << scalarZoneFieldName[nameNo][0] + subdomainName << endl;

        volScalarField fieldValueMode
        (
            IOobject
            (
                scalarZoneFieldName[nameNo][0] + subdomainName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless,
            Field<scalar> (mesh.C().size(), Foam::zero())
        );        


        // assign scalar cell value
        forAll(typesName, typeI)
        {
            word fieldFileName (scalarZoneFieldName[nameNo][typeI]);
            fileName dataFile (dataPath/fieldFileName); 
            IFstream dataStream(dataFile);

            forAll(cellZonesTypeIO[typeI], zoneI)
            {                   
                forAll(mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]], cellI)
                {
                    label cellN (mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]][cellI]);
                    dataStream.read(fieldValueMode[cellN]); 
                }
            } 
        }
        
        // assign boundary value
        forAll(fieldValueMode.boundaryField(), patchI)
        {
            fieldValueMode.boundaryFieldRef().set(patchI, 
                fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
        }

        fieldValueMode.write();          
    }


    // write vector field field
    forAll(vectorZoneFieldName, nameNo)
    {
        Info << "vectorZoneFieldName: " << vectorZoneFieldName[nameNo][0] + subdomainName << endl;

        volVectorField fieldValueMode
        (
            IOobject
            (
                vectorZoneFieldName[nameNo][0] + subdomainName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless,
            Field<vector> (mesh.C().size(), Foam::zero())
        );        


        forAll(typesName, typeI)
        {           
            word fieldFileName (vectorZoneFieldName[nameNo][typeI]);
            fileName dataFile (dataPath/fieldFileName); 
            IFstream dataStream(dataFile);

            forAll(cellZonesTypeIO[typeI], zoneI)
            {
                forAll(mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]], cellI)
                {
                    label cellN (mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]][cellI]);
                    dataStream.read(fieldValueMode[cellN].x()); 
                }

                forAll(mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]], cellI)
                {
                    label cellN (mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]][cellI]);
                    dataStream.read(fieldValueMode[cellN].y()); 
                }

                forAll(mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]], cellI)
                {
                    label cellN (mesh.cellZones()[cellZonesTypeIO[typeI][zoneI]][cellI]);
                    dataStream.read(fieldValueMode[cellN].z()); 
                }
            } 
        }
        
        // assign boundary value
        forAll(fieldValueMode.boundaryField(), patchI)
        {
            fieldValueMode.boundaryFieldRef().set(patchI, 
                fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
        }

        fieldValueMode.write();          
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
