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
    argList::addOption
    (
        "modesFileDir",
        "word",
        "Specify the folder of modes file"
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    fileName dataPath (args.getOrDefault<word>("modesFileDir", mesh.time().path()/"SVD"));

    // List<word> modesNumber ({"0", "1", "2", "3", "4"});
    // // List<word> modesNumber ({"0", "1", "2"});
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
    // List<label> modesNumber (customDict.lookup("modesNumber"));
    label modeNumber (customDict.getLabel("modeNumber"));
    List<word> scalarFieldName (customDict.lookup("scalarFieldName"));
    List<word> vectorFieldName (customDict.lookup("vectorFieldName"));
    // blockType is different cases, case name is insert after fieldName and before modesNumber
    // List<label> blockType (customDict.lookup("blockType"));

    List<label> modesNumber(modeNumber);

    for(label i=0; i<modeNumber; ++i)
        modesNumber[i] = i;

    // write scalar modes field
    forAll(scalarFieldName, nameNo)
    {
        Info << "scalarFieldName: " << scalarFieldName[nameNo] << endl;
           
        forAll(modesNumber, No_)
        {
            fileName modeFieldName(scalarFieldName[nameNo] + "_mode" + name(modesNumber[No_]+1));

            volScalarField fieldValueMode
            (
                IOobject
                (
                    modeFieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimless
            );           

            // asign cell value
            fileName dataFile (dataPath/scalarFieldName[nameNo] + "_mode"); 
            label cellN (0);
            label modeNo = 0;

            Info<< "write the mode field: " << scalarFieldName[nameNo] + "_mode" + name(modesNumber[No_]+1) << endl;

            if(isFile(dataFile))
            {                
                IFstream dataStream(dataFile);
                word dataLine;
                token singleData;   

                while(dataStream.getLine(dataLine) && dataLine != word::null)
                {
                    IStringStream dataString (dataLine);

                    while(modeNo <= modesNumber[No_])
                    {
                        dataString.read(singleData);   
                        ++modeNo;                                
                    }
                    
                    fieldValueMode[cellN] = singleData.scalarToken();

                    ++cellN;
                    modeNo=0;
                }
                            
            }  
            else
            {
                Info << "file: " << dataFile << " is not exist!" << endl;
                // break;
            } 
            
            // assign boundary value
            fileName dataFile (dataPath/scalarFieldName[nameNo] + "_mode"); 
            label cellN (0);
            label modeNo = 0;

            Info<< "write the mode field: " << scalarFieldName[nameNo] + "_mode" + name(modesNumber[No_]+1) << endl;

            if(isFile(dataFile))
            {                
                IFstream dataStream(dataFile);
                word dataLine;
                token singleData;   

                while(dataStream.getLine(dataLine) && dataLine != word::null)
                {
                    IStringStream dataString (dataLine);

                    while(modeNo <= modesNumber[No_])
                    {
                        dataString.read(singleData);   
                        ++modeNo;                                
                    }
                    
                    fieldValueMode[cellN] = singleData.scalarToken();

                    ++cellN;
                    modeNo=0;
                }
                            
            }  
            else
            {
                Info << "file: " << dataFile << " is not exist!" << endl;
                // break;
            }


            fieldValueMode.write();      
        }        

    }


    // write vector modes field
    forAll(vectorFieldName, nameNo)
    {
        Info << "vectorFieldName: " << vectorFieldName[nameNo] << endl;
           
        forAll(modesNumber, No_)
        {
            fileName modeFieldName(vectorFieldName[nameNo] + "_mode" + name(modesNumber[No_]+1));

            RectangularMatrix<scalar> vectorMatrix(mesh.C().size()*3, modeNumber, Foam::Zero);       

            fileName dataFile (dataPath/vectorFieldName[nameNo] + "_mode"); 
            label cellN (0);
            label modeNo = 0;

            Info<< "write the mode field: " << vectorFieldName[nameNo] + "_mode" + name(modesNumber[No_]+1) << endl;

            if(isFile(dataFile))
            {                
                IFstream dataStream(dataFile);
                word dataLine;
                token singleData;   

                while(dataStream.getLine(dataLine) && dataLine != word::null)
                {
                    IStringStream dataString (dataLine);

                    while(modeNo <= modeNumber)
                    {
                        dataString.read(singleData);   
                        vectorMatrix(cellN, modeNo) = singleData.scalarToken();
                        ++modeNo;                                
                    }
                    
                    ++cellN;
                    modeNo=0;
                }
                            
            }  
            else
            {
                Info << "file: " << dataFile << " is not exist!" << endl;
                // break;
            } 

            // Info << "vectorMatrix: " << vectorMatrix.sizes() << endl;
            // Info << vectorMatrix << endl;

            volVectorField fieldValueMode
            (
                IOobject
                (
                    modeFieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimless
            );

            forAll(fieldValueMode, cellI)
            {
                fieldValueMode[cellI] = Vector<scalar> (vectorMatrix(cellI,                     modesNumber[No_]),
                                                        vectorMatrix(cellI + mesh.C().size(),   modesNumber[No_]),
                                                        vectorMatrix(cellI + 2*mesh.C().size(), modesNumber[No_]));
            }

            forAll(fieldValueMode.boundaryField(), patchI)
            {
                fieldValueMode.boundaryFieldRef().set(patchI, 
                    fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
            }

            fieldValueMode.write();      

            // Info << fieldValueMode << endl;
        }        

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
