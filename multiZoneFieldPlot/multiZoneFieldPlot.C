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

    // declaration of field file dir
    fileName dataFile;

    // write scalar field field
    forAll(scalarZoneFieldName, nameNo)
    {
        Info << "scalarZoneFieldName: " << scalarZoneFieldName[nameNo] + subdomainName << endl;

        volScalarField fieldValueMode
        (
            IOobject
            (
                scalarZoneFieldName[nameNo] + subdomainName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless,
            Field<scalar> (mesh.C().size(), Foam::zero())
        );        

        forAll(subdomainZonesList, domainTypeI)
        {
            // file dir of field matrix 
            if(dataPathList.size() ==1)
            {
                dataFile = dataPathList[0]/scalarZoneFieldName[nameNo]; 
            }
            else
            {
                dataFile = dataPathList[domainTypeI]/scalarZoneFieldName[nameNo]; 
            }
            
            if(isFile(dataFile))
            {                
                IFstream dataStream(dataFile);

                // output the subdomains List
                IOList<List<label>> subdomainsIO
                (
                    IOobject
                    (
                        subdomainZonesList[domainTypeI],
                        runTime.caseConstant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
                
                // asign cell value
                forAll(subdomainsIO, subdomainI)
                {  
                    forAll(subdomainsIO[subdomainI], zoneI)
                    {                   
                        forAll(mesh.cellZones()[subdomainsIO[subdomainI][zoneI]], cellI)
                        {
                            label cellN (mesh.cellZones()[subdomainsIO[subdomainI][zoneI]][cellI]);
                            dataStream.read(fieldValueMode[cellN]); 
                        }
                    }
                }   
            }  
            else
            {
                Info << "file: " << dataFile << " is not exist!" << endl;
                // break;
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
        Info << "vectorZoneFieldName: " << vectorZoneFieldName[nameNo] + subdomainName << endl;

        volVectorField fieldValueMode
        (
            IOobject
            (
                vectorZoneFieldName[nameNo] + subdomainName,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless,
            Field<vector> (mesh.C().size(), Foam::zero())
        );        


        forAll(subdomainZonesList, domainTypeI)
        {
            // file dir of field matrix 
            if(dataPathList.size() ==1)
            {
                dataFile = dataPathList[0]/vectorZoneFieldName[nameNo]; 
            }
            else
            {
                dataFile = dataPathList[domainTypeI]/vectorZoneFieldName[nameNo]; 
            }

            if(isFile(dataFile))
            {                
                IFstream dataStream(dataFile);

                // output the subdomains List
                IOList<List<label>> subdomainsIO
                (
                    IOobject
                    (
                        subdomainZonesList[domainTypeI],
                        runTime.caseConstant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                // asign cell value
                forAll(subdomainsIO, subdomainI)
                {  
                    forAll(subdomainsIO[subdomainI], zoneI)
                    {                   
                        forAll(mesh.cellZones()[subdomainsIO[subdomainI][zoneI]], cellI)
                        {
                            label cellN (mesh.cellZones()[subdomainsIO[subdomainI][zoneI]][cellI]);
                            dataStream.read(fieldValueMode[cellN].x()); 
                        }
                    }

                    forAll(subdomainsIO[subdomainI], zoneI)
                    {                   
                        forAll(mesh.cellZones()[subdomainsIO[subdomainI][zoneI]], cellI)
                        {
                            label cellN (mesh.cellZones()[subdomainsIO[subdomainI][zoneI]][cellI]);
                            dataStream.read(fieldValueMode[cellN].y()); 
                        }
                    }

                    forAll(subdomainsIO[subdomainI], zoneI)
                    {                   
                        forAll(mesh.cellZones()[subdomainsIO[subdomainI][zoneI]], cellI)
                        {
                            label cellN (mesh.cellZones()[subdomainsIO[subdomainI][zoneI]][cellI]);
                            dataStream.read(fieldValueMode[cellN].z()); 
                        }
                    }
                }                 
            }  
            else
            {
                Info << "file: " << dataFile << " is not exist!" << endl;
                // break;
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
