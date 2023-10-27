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

// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone

#include "fvCFD.H"
#include "IFstream.H"
#include "stringOps.H"

int main(int argc, char *argv[])
{
    // add timeSelector to generate all ROM results
    timeSelector::addOptions();

    argList::addOption
    (
        "modesFileDir",
        "word",
        "Specify the folder of modes file"
    );

    // Rsh, add dict select function
    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative modePlotDict"
    );
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    // write the field value to specific time folder
    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs.last(), 0); 

    fileName dataPath (args.getOrDefault<word>("modesFileDir", runTime.globalPath()/"SVD"));

    // Rsh, the default dict name or input selector
    //  - must plus the dict select option in the beginning of the code
    const word dictName("DGReconDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;
    IOdictionary DGReconDict(dictIO);

    // read data from dictionary
    // List<label> modesNumber (customDict.lookup("modesNumber"));
    label elementNum (DGReconDict.getLabel("elementNum"));
    List<word> scalarFileName (DGReconDict.lookup("scalarFileName"));
    List<word> vectorFileName (DGReconDict.lookup("vectorFileName"));
    // blockType is different cases, case name is insert after fieldName and before modesNumber
    // List<label> blockType (customDict.lookup("blockType"));

    RectangularMatrix<scalar> scalarFieldMat(mesh.cellZones()[0].size(), elementNum);
    RectangularMatrix<scalar> vectorFieldMat(3*mesh.cellZones()[0].size(), elementNum);

    // write scalar field
    forAll(scalarFileName, nameNo)
    {
        fileName dataFile (dataPath/scalarFileName[nameNo]);
        if(isFile(dataFile))
        {                
            IFstream dataStream(dataFile);
            word dataLine;
            label row(0);

            while(dataStream.getLine(dataLine) && dataLine != word::null)
            {
                IStringStream dataString (dataLine);
                token singleData;  // token stores the data read from IFstream 

                for(label elementI = 0; elementI < elementNum; ++elementI)
                {
                    dataString.read(singleData);    
                    scalarFieldMat(row, elementI) = singleData.scalarToken();
                }   
                ++row;
            }                       
        }  
        else
        {
            Info << "file: " << dataFile << " is not exist!" << endl;
            // break;
        }

        // create mode field by copying T
        volScalarField sField
        (
            IOobject
            (
                scalarFileName[nameNo],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
        );

        // assign cell value
        forAll(mesh.cellZones(), ZoneI)
        {       
            forAll(mesh.cellZones()[ZoneI], cellI)
            {       
                sField[mesh.cellZones()[ZoneI][cellI]] = scalarFieldMat(cellI, ZoneI);
            }
        }

        // assign boundary value
        forAll(sField.boundaryField(), patchI)
        {
            sField.boundaryFieldRef().set(patchI, 
                fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], sField));
        }

        sField.write(); 
    }

    // write vector modes field
    forAll(vectorFileName, nameNo)
    {
        fileName dataFile (dataPath/vectorFileName[nameNo]);
        if(isFile(dataFile))
        {                
            IFstream dataStream(dataFile);
            word dataLine;
            label row(0);

            while(dataStream.getLine(dataLine) && dataLine != word::null)
            {
                IStringStream dataString (dataLine);
                token singleData;  // token stores the data read from IFstream 

                for(label elementI = 0; elementI < elementNum; ++elementI)
                {
                    dataString.read(singleData);    
                    vectorFieldMat(row, elementI) = singleData.scalarToken();
                }   
                ++row;
            }                       
        }  
        else
        {
            Info << "file: " << dataFile << " is not exist!" << endl;
            // break;
        }

        volVectorField vField
        (
            IOobject
            (
                vectorFileName[nameNo],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
        );

        // assign cell value
        forAll(mesh.cellZones(), ZoneI)
        {       
            forAll(mesh.cellZones()[ZoneI], cellI)
            {
                vField[mesh.cellZones()[ZoneI][cellI]] = Vector<scalar>(vectorFieldMat(cellI,                     ZoneI),
                                                                        vectorFieldMat(cellI + mesh.cellZones()[ZoneI].size(),   ZoneI),
                                                                        vectorFieldMat(cellI + 2*mesh.cellZones()[ZoneI].size(), ZoneI));
            }
        }

        forAll(vField.boundaryField(), patchI)
        {
            vField.boundaryFieldRef().set(patchI, 
                fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], vField));
        }

        vField.write();      
    }        


    Info << "\nEnd\n" << endl;

    return 0;

}


// ************************************************************************* //
