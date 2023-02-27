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

// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone

#include "fvCFD.H"
#include "wordRe.H"

int main(int argc, char *argv[])
{
    // Rsh, add dict select function
    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative cellZoneDict"
    );

    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // add dict Selector option
    const word dictName("cellZoneDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

    // create dictionary object and read the dictionary
    IOdictionary customDict(dictIO);
    #include "readDict.H"


    List<word> allBoundaryNames(customDict.lookup("allBoundaryNames"));


    word outputDir(customDict.getWord("outputDir"));
    // create the SVD folder if it is not exist.
    if(!isDir(runTime.globalPath()/outputDir))
        mkDir(runTime.globalPath()/outputDir);

    fileName dataFile;


    // read the fields
    List<word> scalarFields(customDict.lookup("scalarFields"));
    List<word> vectorFields(customDict.lookup("vectorFields"));


    fileName snapshotFileName;
    autoPtr<OFstream> outputFilePtr;

    // scalarFields cell value
    forAll(scalarFields, fieldI)
    {   
        forAll(typesName, typeI)
        {
            label totalzones(0);
            forAll(cellZonesType[typeI], caseI)
            {
                totalzones += cellZonesType[typeI][caseI].size(); 
            }

            label zoneCells(returnReduce(mesh.cellZones()[cellZonesType[typeI][0][0]].size(), sumOp<label>()) );
            Info<< "Size of the snapshots, " << scalarFields[fieldI]+name(typeI) << ", " 
                << zoneCells << ", " << totalzones << endl;

            snapshotFileName = scalarFields[fieldI]+name(typeI);
            if (Pstream::master())
                outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

            for(label caseI=0; caseI<casesNumber; ++caseI)
            {
                Info<< "read field value: " << scalarFields[fieldI]+name(typeI)+"_"+name(CasesID[caseI]) << endl;
                volScalarField tmpField
                (
                    IOobject
                    (
                        scalarFields[fieldI]+"_"+name(CasesID[caseI]),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                forAll(cellZonesType[typeI][caseI], zoneI)
                {                   
                    // collected cell value from individual processors
                    List<scalar> cellValue;
                    forAll(mesh.cellZones()[cellZonesType[typeI][caseI][zoneI]], cellI)
                    {
                        label cellN (mesh.cellZones()[cellZonesType[typeI][caseI][zoneI]][cellI]);
                        cellValue.append(tmpField[cellN]);
                    }
                    List< List<scalar> > gatheredCellValue(Pstream::nProcs());
                    gatheredCellValue[Pstream::myProcNo()] = cellValue;
                    Pstream::gatherList(gatheredCellValue);
                    List<scalar> cellValue_  = ListListOps::combine<List<scalar> >(gatheredCellValue, accessOp<List<scalar> >());

                    // Loop the collected cellValue list, "write" to file can only run in the master processor to avoid error
                    // -- while "Info" do not need to run in the master processor
                    if (Pstream::master())
                    {
                        forAll(cellValue_, cellI)
                        {
                            outputFilePtr().write(cellValue_[cellI]) << " ";
                        }
                        outputFilePtr() << endl;  
                    }
                } 
            }
        }
    }


    // vectorFields cell value
    forAll(vectorFields, fieldI)
    {   
        forAll(typesName, typeI)
        {           
            label totalzones(0);
            forAll(cellZonesType[typeI], caseI)
            {
                totalzones += cellZonesType[typeI][caseI].size(); 
            }

            label zoneCells(returnReduce(mesh.cellZones()[cellZonesType[typeI][0][0]].size(), sumOp<label>()) );
            Info<< "Size of the snapshots, " << vectorFields[fieldI]+name(typeI) << ", "
                << zoneCells << ", " << totalzones << endl;

            snapshotFileName = vectorFields[fieldI]+name(typeI);
            if (Pstream::master())
                outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

            
            for(label caseI=0; caseI<casesNumber; ++caseI)
            {
                Info<< "read field value: " << vectorFields[fieldI]+name(typeI)+"_"+name(CasesID[caseI]) << endl;
                volVectorField tmpField
                (
                    IOobject
                    (
                        vectorFields[fieldI]+"_"+name(CasesID[caseI]),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                forAll(cellZonesType[typeI][caseI], zoneI)
                {
                    // collected cell value from individual processors
                    List<vector> cellValue;
                    forAll(mesh.cellZones()[cellZonesType[typeI][caseI][zoneI]], cellI)
                    {
                        label cellN (mesh.cellZones()[cellZonesType[typeI][caseI][zoneI]][cellI]);
                        cellValue.append(tmpField[cellN]);
                    }
                    List< List<vector> > gatheredCellValue(Pstream::nProcs());
                    gatheredCellValue[Pstream::myProcNo()] = cellValue;
                    Pstream::gatherList(gatheredCellValue);
                    List<vector> cellValue_  = ListListOps::combine<List<vector> >(gatheredCellValue, accessOp<List<vector> >());

                    // Loop the collected cellID list
                    if (Pstream::master())
                    {                      
                        forAll(cellValue_, cellI)
                        {
                            outputFilePtr() << cellValue_[cellI].x() << " ";
                        }
                        forAll(cellValue_, cellI)
                        {
                            outputFilePtr() << cellValue_[cellI].y() << " ";
                        }
                        forAll(cellValue_, cellI)
                        {
                            outputFilePtr() << cellValue_[cellI].z() << " ";
                        }
                        outputFilePtr() << endl;    
                    }
                } 
            }
        }
    }
    

    Info << "\nEnd\n" << endl;

    return 0;

}


// ************************************************************************* //
