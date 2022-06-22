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
// #include "emptyPolyPatch.H"
#include "IFstream.H"
#include "stringOps.H"
#include "Time.H"

int main(int argc, char *argv[])
{
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    fileName dataPath (mesh.time().path()/"postProcessing");

    // const word dictName("reconstructDict");
    const word dictName("svdDict");
    // Create and input-output object - this holds the path to the dict and its name
    IOdictionary reconstructDict
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
    // timeValueNum, total time interval number 
    // reconstructTime, time point need output 
    // reconstructFileName, the file name of ROM result
    label timeValueNum (reconstructDict.lookupOrDefault<label>("timeValueNum", 10));
    List<label> reconstructTime (reconstructDict.lookup("reconstructTime"));
    word reconstructFileName (mesh.time().path()/"SVD"/"snapshots_calculate");
    scalar writeTime_ (reconstructDict.lookupOrDefault<label>("writeTime", runTime.endTime().value()));

    // The reconstruct Matrix
    RectangularMatrix<scalar> romResultMatrix(mesh.C().size(), timeValueNum);
    Info<< "romResultMatrix: " << romResultMatrix.sizes() << endl;
    word romResultMatrixName (mesh.time().path()/"SVD"/"romResultMatrixName");
    OFstream romResultMatrixOutput (romResultMatrixName);
    label row_ = 0;
    label column_ = 0;

    // read ROM result
    IFstream resultStream(reconstructFileName);
    word rowValue;


    while(row_ < mesh.C().size() && ! resultStream.eof())
    {
        resultStream.getLine(rowValue);        
        IStringStream resultLineStream(rowValue);
        
        // Info<< "test " << endl;
        // Info<< rowValue << endl;

        column_ = 0;
        while (column_ < timeValueNum && ! resultLineStream.eof())
        {
            resultLineStream.read(romResultMatrix[row_][column_]);

            romResultMatrixOutput.width(16);
            romResultMatrixOutput << romResultMatrix[row_][column_];

            ++column_;
            // Info<< column_ << endl;            
        }
        romResultMatrixOutput << endl;
        ++row_;

        // Info<< row_ << endl;
    }

    // Info<< "row: " << row_ << endl 
    //     << "column: " << column_ << endl;
    
    // set the time folder to write the field value
    // scalar writeTime_(40);
    runTime.setTime(instant(writeTime_), 5);

    forAll(reconstructTime, timeI)
    {
        volScalarField fieldValue
        (
            IOobject
            (
                "reconstruct_" + name(reconstructTime[timeI]),
                mesh.time().path(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
        );

        forAll(fieldValue, cellI)
        {
            fieldValue[cellI] = romResultMatrix[cellI][reconstructTime[timeI]];
        }

        fieldValue.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
