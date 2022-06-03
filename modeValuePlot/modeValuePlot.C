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
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    fileName dataPath (mesh.time().path()/"postProcessing");

    // List<word> modeNumber ({"0", "1", "2", "3", "4"});
    // // List<word> modeNumber ({"0", "1", "2"});
    // List<word> fieldName ({"Ux", "Uy", "Uz", "magU", "gradpx", "gradpy", "gradpz"});

    const word dictName("modePlotDict");

    // Create and input-output object - this holds the path to the dict and its name
    IOdictionary customDict
    (
        IOobject
        (
            dictName, // name of the file
            mesh.time().system(), // path to where the file is
            mesh, // reference to the mesh needed by the constructor
            IOobject::MUST_READ // indicate that reading this dictionary is compulsory
        )
    );

    // read data from dictionary
    List<label> modeNumber (customDict.lookup("modeNumber"));
    List<word> fieldName (customDict.lookup("fieldName"));
    List<label> blockType (customDict.lookup("blockType"));

    forAll(fieldName, nameNo)
    {
        Info << "FieldName: " << fieldName[nameNo] << endl;

        // // read zone name of the field field value
        // IFstream zoneNameStream(dataPath/fieldName[nameNo]);
        // List<word> zoneName (mesh.cellZones().size());

        // // use IStringStream to read zoneName and put them into a List
        // label zoneNumber (0);
        // word zoneNameLine;
        // zoneNameStream.getLine(zoneNameLine);
        // IStringStream zoneNameString (zoneNameLine);

        // while (! zoneNameString.eof())
        // {
        //     zoneNameString >> zoneName[zoneNumber];
        //     // Info << "zoneName: " << zoneName[zoneNumber] << endl;
        //     zoneNumber += 1;
        // }
        // zoneName.resize(zoneNumber);
           
        forAll(modeNumber, No_)
        {
            fileName modeFieldName(fieldName[nameNo] + "_mode" + name(modeNumber[No_]));

            Field<scalar> zeroScalarField (mesh.C().size(), Foam::zero());
            volScalarField fieldValueMode
            (
                IOobject
                (
                    modeFieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimVelocity,
                zeroScalarField
            );

            forAll(blockType, blockTypeNo)
            {
                // read zone name of the field field value
                word blockFieldName = fieldName[nameNo] + "_" + name(blockType[blockTypeNo]);
                Info << "blockFieldName: " << blockFieldName << endl;
                IFstream zoneNameStream(dataPath/blockFieldName);
                List<word> zoneName (mesh.cellZones().size());

                // use IStringStream to read zoneName and put them into a List
                label zoneNumber (0);
                word zoneNameLine;
                zoneNameStream.getLine(zoneNameLine);
                // Info << "zoneNameLine: " << zoneNameLine << endl;
                IStringStream zoneNameString (zoneNameLine);

                while (! zoneNameString.eof())
                {
                    zoneNameString >> zoneName[zoneNumber];
                    // Info << "zoneName: " << zoneName[zoneNumber] << endl;
                    zoneNumber += 1;
                }
                zoneName.resize(zoneNumber);
                
                // Info << "test1 " << endl;

                fileName dataFile (dataPath/fieldName[nameNo]+ 
                                    "_" + name(blockType[blockTypeNo]) + 
                                    "_mode" + name(modeNumber[No_]));                

                // Info << "dataFile: " << dataFile << endl
                //      << "modeFieldName: " << modeFieldName << endl;

                // Field<scalar> zeroScalarField (mesh.C().size(), Foam::zero());

                if(isFile(dataFile))
                {
                    // volScalarField fieldValueMode
                    // (
                    //     IOobject
                    //     (
                    //         modeFieldName,
                    //         mesh.time().timeName(),
                    //         mesh,
                    //         IOobject::READ_IF_PRESENT,
                    //         IOobject::AUTO_WRITE
                    //     ),
                    //     mesh,
                    //     dimVelocity,
                    //     zeroScalarField
                    // );
                    
                    IFstream dataStream(dataFile);

                    label firstZoneI = mesh.cellZones().findZoneID(zoneName[0]);

                    // Info << "firstZoneI: " << firstZoneI << endl;

                    forAll(mesh.cellZones()[firstZoneI], cellI)
                    {
                        forAll(zoneName, zoneNameI)
                        {
                            label zoneI = mesh.cellZones().findZoneID(zoneName[zoneNameI]);
                            label cell = mesh.cellZones()[zoneI][cellI];
                            dataStream.read(fieldValueMode[cell]);            
                        }
                    }

                    fieldValueMode.write();
                }  
            }   
            
        }

        

        

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
