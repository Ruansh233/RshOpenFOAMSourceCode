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
#include "IFstream.H"
#include "stringOps.H"

int main(int argc, char *argv[])
{
    // add timeSelector to generate all ROM results
    timeSelector::addOptions();
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    fileName dataPath (mesh.time().path()/"SVD");
    autoPtr<OFstream> outputFilePtr;

    const word dictName("svdDict");

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
    label snapshotsNum (customDict.getLabel("snapshotsNum"));
    label modesNum (customDict.getLabel("modeExtract"));
    List<word> calCoeffFile (customDict.lookup("calCoeffFile"));
    List<word> modesMatrixName (customDict.lookup("modesName"));
    List<word> snapshotsMatrixName (customDict.lookup("snapshotsMatrixName"));

    // variable to count the structure of temporal Coefficient
    label snapshotsN(0);
    label modesN(0);

    // // read calculated temporal Coefficient and store into a matrix 
    // forAll(calCoeff, subdomainI)
    // {
        // calculated temporal Coefficient file name
        // fileName dataFile (dataPath/calCoeff[subdomainI]); 
        fileName dataFile (dataPath/calCoeffFile[0]);
        // the matrix store Coefficient
        RectangularMatrix<scalar> calModecoeff(snapshotsNum, modesNum);            

        if(isFile(dataFile))
        {                
            IFstream dataStream(dataFile);
            word dataLine;

            while(snapshotsN < snapshotsNum)
            {
                dataStream.getLine(dataLine);
                IStringStream dataString (dataLine);
                token singleData;   

                while(!dataString.eof())
                {
                    dataString.read(singleData);    
                    if(singleData.isScalar())
                    {
                        calModecoeff(snapshotsN, modesN) = singleData.scalarToken();
                        ++modesN;                            
                    }         
                }
                ++snapshotsN;
                modesN = 0;
            }                          
        }  
        else
        {
            Info << "file: " << dataFile << " is not exist!" << endl;
            // break;
        }

        Info<< "calculated Coefficient matrix is: " << calModecoeff.m() << " * " 
        << calModecoeff.n() << endl; 

        if(snapshotsN != snapshotsNum)
            Info<< "The rows of calModecoeff should equal to snapshots numbers" << endl;

        // // write calModecoeff
        // dataFile = mesh.time().path()/"SVD"/"calModecoeff_test";    
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < calModecoeff.m(); ++row)
        // {
        //     for (label column = 0; column < calModecoeff.n(); ++column)
        //     {
        //         outputFilePtr().width(12);
        //         outputFilePtr() << calModecoeff(row, column);
        //     }
        //     outputFilePtr() << endl;
        // }


        // read real snapshots matrix
        fileName snapshotsMFile (dataPath/snapshotsMatrixName[0]);
        RectangularMatrix<scalar> snapshotsM(mesh.C().size(), snapshotsNum);            

        label cellN (0);
        snapshotsN = 0;

        if(isFile(snapshotsMFile))
        {                
            IFstream dataStream(snapshotsMFile);
            word dataLine;

            while(dataStream.getLine(dataLine) && dataLine != word::null)
            {
                IStringStream dataString (dataLine);
                token singleData;   

                while(!dataString.eof())
                {
                    dataString.read(singleData);
                    if(singleData.isScalar())
                    {
                        snapshotsM(cellN, modesN) = singleData.scalarToken();
                        ++modesN;                            
                    }     
                    if(singleData.isLabel())
                    {
                        snapshotsM(cellN, modesN) = scalar(singleData.labelToken());
                        ++modesN;                   
                    }                    
                }
                ++cellN;
                modesN = 0;
            }
                          
        }  
        else
        {
            Info << "file: " << snapshotsMFile << " is not exist!" << endl;
            // break;
        }
        Info<< "snapshots matrix is: " << snapshotsM.m() << " * " 
            << snapshotsM.n() << endl; 

        // // write snapshotsM
        // dataFile = mesh.time().path()/"SVD"/"snapshotsM_test";    
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < snapshotsM.m(); ++row)
        // {
        //     for (label column = 0; column < snapshotsM.n(); ++column)
        //     {
        //         outputFilePtr().width(12);
        //         outputFilePtr() << snapshotsM(row, column);
        //     }
        //     outputFilePtr() << endl;
        // }

        // read mode matrix
        // only the required number of modes are read
        fileName modeMatrixFile (dataPath/modesMatrixName[0]);
        RectangularMatrix<scalar> modesMatrix(mesh.C().size(), modesNum);            

        cellN = 0;
        modesN = 0;

        if(isFile(modeMatrixFile))
        {                
            IFstream dataStream(modeMatrixFile);
            word dataLine;

            while(dataStream.getLine(dataLine) && dataLine != word::null)
            {
                IStringStream dataString (dataLine);
                token singleData;   

                while(modesN < modesNum)
                {
                    dataString.read(singleData);
                    if(singleData.isScalar())
                    {
                        modesMatrix(cellN, modesN) = singleData.scalarToken();
                        ++modesN;                            
                    }                    
                }

                ++cellN;
                modesN = 0;                
            }
                          
        }  
        else
        {
            Info << "file: " << modeMatrixFile << " is not exist!" << endl;
            // break;
        }
        Info<< "modes matrix is: " << modesMatrix.m() << " * " 
            << modesMatrix.n() << endl; 
        
        // // write modesMatrix
        // dataFile = mesh.time().path()/"SVD"/"modesMatrix_test";    
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < modesMatrix.m(); ++row)
        // {
        //     for (label column = 0; column < modesMatrix.n(); ++column)
        //     {
        //         outputFilePtr().width(12);
        //         outputFilePtr() << modesMatrix(row, column);
        //     }
        //     outputFilePtr() << endl;
        // }


        // snapshots of ROM and the error estimate 
        RectangularMatrix<scalar> calSnapshotsM(mesh.C().size(), snapshotsNum); 
        calSnapshotsM = modesMatrix * calModecoeff.T();

        Info<< "calculated snapshots matrix is: " << calSnapshotsM.m() << " * " 
            << calSnapshotsM.n() << endl; 

        // // write calSnapshotsM
        // dataFile = mesh.time().path()/"SVD"/"calSnapshotsM_test";    
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < calSnapshotsM.m(); ++row)
        // {
        //     for (label column = 0; column < calSnapshotsM.n(); ++column)
        //     {
        //         outputFilePtr().width(12);
        //         outputFilePtr() << calSnapshotsM(row, column);
        //     }
        //     outputFilePtr() << endl;
        // }

        RectangularMatrix<scalar> errorMatrix(mesh.C().size(), snapshotsNum); 
        errorMatrix = calSnapshotsM - snapshotsM;

        // // write calModecoeff
        // dataFile = mesh.time().path()/"SVD"/"errorMatrix_test";    
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < errorMatrix.m(); ++row)
        // {
        //     for (label column = 0; column < errorMatrix.n(); ++column)
        //     {
        //         outputFilePtr().width(12);
        //         outputFilePtr() << errorMatrix(row, column);
        //     }
        //     outputFilePtr() << endl;
        // }
        
        List<scalar> averageError(snapshotsNum);
        List<scalar> maxError(snapshotsNum);
        List<scalar> tempValue(mesh.C().size());

        for(label columns = 0; columns < errorMatrix.n(); ++columns)
        {
            for(label rows = 0; rows < errorMatrix.m(); ++rows)
            {
                tempValue[rows] = mag(errorMatrix(rows, columns)/snapshotsM(rows, columns));
            }
            maxError[columns] = max(tempValue);        
            averageError[columns] = average(tempValue);
        } 


        // output the error
        dataFile = mesh.time().path()/"SVD"/"errorList";    
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < averageError.size(); ++row)
        {
            outputFilePtr().width(16);
            outputFilePtr() << row;
            outputFilePtr().width(16);
            outputFilePtr() << averageError[row];
            outputFilePtr().width(16);
            outputFilePtr() << maxError[row];

            outputFilePtr() << endl;
        }

    // write data of specific time
    // read time list
    // List<label> outputTime (customDict.lookup("outputTime"));
    instantList timeDirs = timeSelector::select0(runTime, args);
    scalar timeInterval (customDict.getScalar("timeInterval"));

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], 1);        
        Info<< "runtime: " << runTime.timeName() << endl;

        volScalarField tempFieldValue
        (
            IOobject
            (
                "T",
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        volScalarField calFieldValue
        (
            IOobject
            (
                "calFieldValue",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            tempFieldValue
        );

        label timeColumn (runTime.value()/timeInterval);
        if(timeColumn >= calSnapshotsM.n())
            --timeColumn;

        forAll(calFieldValue, cellI)
        {
            calFieldValue[cellI] = calSnapshotsM(cellI, timeColumn);
        }

        calFieldValue.write();
    }

    // }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
