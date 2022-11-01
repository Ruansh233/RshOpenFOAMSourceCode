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

    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    if (mesh.cellZones().size() == 0)
    {
        FatalErrorIn("cellZoneFieldValueROM")
            << "There is no cellZone in this mesh"
            << exit(FatalError);
    }

    fileName dataPath (mesh.time().path()/"SVD");

    label elementNum (5);

    // read the matrix
    // read cell modes matrix
    fileName dataFile1 (dataPath/"calSnapshotsM");
    fileName dataFile2 (dataPath/"errorM");

    RectangularMatrix<scalar> calSnapshotsM(mesh.C().size(), elementNum);
    RectangularMatrix<scalar> errorM(mesh.C().size(), elementNum);

    if(isFile(dataFile1))
    {                
        IFstream dataStream(dataFile1);
        word dataLine;
        label row(0);

        while(dataStream.getLine(dataLine) && dataLine != word::null)
        {
            IStringStream dataString (dataLine);
            token singleData;  // token stores the data read from IFstream 

            for(label elementI = 0; elementI < elementNum; ++elementI)
            {
                dataString.read(singleData);    
                calSnapshotsM(row, elementI) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile1 << " is not exist!" << endl;
        // break;
    }

    if(isFile(dataFile2))
    {                
        IFstream dataStream(dataFile2);
        word dataLine;
        label row(0);

        while(dataStream.getLine(dataLine) && dataLine != word::null)
        {
            IStringStream dataString (dataLine);
            token singleData;  // token stores the data read from IFstream 

            for(label elementI = 0; elementI < elementNum; ++elementI)
            {
                dataString.read(singleData);    
                errorM(row, elementI) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile2 << " is not exist!" << endl;
        // break;
    }

    // ---------------------------------------------------------------------------------//
    // assign cell value to each cellZone
    // create field value
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    // create mode field by copying T
    volScalarField calT
    (
        IOobject
        (
            "calT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T
    );

    // create mode field by copying T
    volScalarField error
    (
        IOobject
        (
            "error",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T
    );


    // assign cell value
    forAll(mesh.cellZones(), ZoneI)
    {       
        forAll(mesh.cellZones()[ZoneI], cellI)
        {
             
            calT[mesh.cellZones()[ZoneI][cellI]] = calSnapshotsM(cellI, ZoneI);
            error[mesh.cellZones()[ZoneI][cellI]] = errorM(cellI, ZoneI);
        }
    }

    
    // write the field value to specific time folder
    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs.last(), 0); 

    calT.write();
    error.write();

    Info << "\nEnd\n" << endl;

    return 0;

}


// ************************************************************************* //
