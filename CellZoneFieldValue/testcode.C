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
#include "cellSet.H"
// #include "objectRegistry.H"

// #include "readFieldTemplate.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // Get access to a custom dictionary
    // dictionary customDict;
    const word dictName("cellZoneDict");

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

    // customDict = IOdictionary(dictIO);

    word runTimeName_ (customDict.lookup("readTime"));
    word fieldName_ (customDict.lookup("fieldName"));
    // Info<< "field time is: " << runTimeName_ << endl << endl;

    // word cellSetName_ (customDict.lookup("cellSetName"));   
    
    // List<word> cellSetsName_ (customDict.lookup("cellSetName"));
    List<word> cellZonesName_ (customDict.lookup("cellZoneName"));


    // if (fieldName_ == "U")
    // {
    //     Info<< "field time is: " << fieldName_ << endl << endl;
    // }


    // Info << "mesh.time().path(): " << mesh.time().path() << endl;
    // Info << "runTime.system(): " << runTime.system() << endl;
    // // Info << "args: " << args << endl;
    // Info << "args.rootPath(): " << args.rootPath() << endl;
    // Info << "args.caseName(): " << args.caseName() << endl;

    // mesh.time().path(): "/mnt/d/Learning/OpenFOAMsourcecode/CellZoneFieldValue/testCase"
    // runTime.system(): system
    // args.rootPath(): "/mnt/d/Learning/OpenFOAMsourcecode/CellZoneFieldValue"
    // args.caseName(): "testCase"

    // word myName (customDict.subDict("hello").lookup("name"));
    // Info << "maName is: " << myName << endl;


    // readField ( runTimeName_, fieldName_, mesh );

    
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField p // note that pressure is a scalar field
	(
		IOobject
		(
		    "p", // name of the field
		    runTimeName_, // name of the current time, i.e. the time folder to read from
		    mesh,
		    IOobject::MUST_READ, // always gets imported, will throw an error if the field is missing
		    IOobject::NO_WRITE // will get saved automatically when the controlDict parameters will request it
		),
		mesh // initialises the field to match the size of the mesh with default (0) values
	);

    volVectorField gradp(fvc::grad(p));

    // volVectorField gradU(fvc::grad(U.component(0)));

    // volScalarField Ux_(U.component(vector::X));
    // volVectorField gradUx_(fvc::grad(U.component(vector::X)));

    // Info << "X: " << vector::X << endl;

    // Info << fvc::grad;

    // volScalarField magU(mag(U));
    // volVectorField gradU(fvc::grad(U));

    // volScalarField p // note that pressure is a scalar field
	// (
	// 	IOobject
	// 	(
	// 	    "p", // name of the field
	// 	    runTimeName_, // name of the current time, i.e. the time folder to read from
	// 	    mesh,
	// 	    IOobject::MUST_READ, // always gets imported, will throw an error if the field is missing
	// 	    IOobject::NO_WRITE // will get saved automatically when the controlDict parameters will request it
	// 	),
	// 	mesh // initialises the field to match the size of the mesh with default (0) values
	// );

    // volScalarField magU // note that pressure is a scalar field
	// (
	// 	IOobject
	// 	(
	// 	    "mag(U)", // name of the field
	// 	    runTimeName_, // name of the current time, i.e. the time folder to read from
	// 	    mesh,
	// 	    IOobject::MUST_READ, // always gets imported, will throw an error if the field is missing
	// 	    IOobject::NO_WRITE // will get saved automatically when the controlDict parameters will request it
	// 	),
	// 	mesh // initialises the field to match the size of the mesh with default (0) values
	// );

    List<scalarList> outputData_(cellZonesName_.size());
    // List<scalarList> outputData_(cellSetsName_.size());


    // scalarList testList_;
    // testList_.append(5);

    // Info << testList_[0] << endl;

    // outputData_[0][0] = 5.0;

    forAll(cellSetsName_, setID)
    {
        
        Info << endl << "cellSet name: " << cellSetsName_[setID] << endl;

        // Info << endl << cellSetsName_.size() << endl;
        
        word cellSetName_ (cellSetsName_[setID]);

        cellSet selectedCells(mesh, cellSetName_);   // define a cellset based on mesh and cellset name
        const labelList& cellSets_ = selectedCells.toc(); //get cell ID in a cellset

        // forAll(cellSets_, i)
        // {
        //     const label cellI = cellSets_[i];
        //     Info << "The center of cell " << cellI << " is " << mesh.C()[cellI]
        //          << ", velocity is " << gradUx_[cellI] << ", velocity is " << U[cellI] << endl;
        // }

        // outputData_[setID] = U.component(setID)

        forAll(cellSets_, i)
        {
            const label cellI = cellSets_[i];
            outputData_[setID].append(U[cellI].z());

            // Info << "outputData_[setID][cellI]: " << outputData_[setID][cellI] << endl;
            // outputFile_ << U[cellI].x();
            // Info << mesh.C()[cellI].x() << " " << mesh.C()[cellI].y() << " " << mesh.C()[cellI].z() << endl;
            // Info << "The center of cell " << cellI << " is " << mesh.C()[cellI] << endl
            //      << ", pressure gradient is " << gradp[cellI] << ", pressure is " << p[cellI] << endl;
        }
    }

    // Info << outputData_.size() << endl;
    // Info << outputData_[0].size() << endl;

    // Create the output path directory
    fileName outputDir = mesh.time().path()/"postProcessing";
    // Info << mesh.time().path() << nl;
    // Creathe the directory
    mkDir(outputDir);

    // // File pointer to direct the output to
    // autoPtr<OFstream> outputFilePtr;
    // // Open the file in the newly created directory
    // outputFilePtr.reset(new OFstream(outputDir/"customOutputFile.dat"));

    OFstream outputFile_(outputDir/fieldName_);

    forAll(outputData_[0], dataID)
    {
        forAll(outputData_, blockID)
        {
            outputFile_.width(16);
            outputFile_ << outputData_[blockID][dataID];
        }
        outputFile_ << endl;
    }

    

    // Write stuff
    // outputFilePtr() << "# This is a header" << endl;
    // outputFilePtr() << "0 1 2 3 4 5" << endl;

    // Info << outputData_ << endl;

    // label testLabel_;
    // testLabel_ = 5;

    // Info << endl << testLabel_ << endl;

    // forAll( p.boundaryField(), ipatch ) 
    // {
    //     Info << endl << "Patch " << ipatch << ": " << mesh.boundary()[ipatch].name() << endl
    //          << "patch type is: " << p.boundaryField()[ipatch];
    //     /// each face in the the patch
    //     forAll( p.boundaryField()[ipatch], iface ) 
    //     {
    //         Info << p.boundaryField()[ipatch][iface] << endl;
    //         /// ...
    //     }
    // }

    // Rsh, 
    // U[cellI].x() is the x component of the velocity vector

    // word cellSetName_ (cellSetsName_[0]);

    // cellSet selectedCells(mesh, cellSetName_);
    // const labelList& cellSets_ = selectedCells.toc();

    
    // Info<< "field time is: " << runTimeName_ << endl << endl;

    // // Info<< "Reading field p\n" << endl;
	// // volScalarField p // note that pressure is a scalar field
	// // (
	// // 	IOobject
	// // 	(
	// // 	    "p", // name of the field
	// // 	    runTimeName_, // name of the current time, i.e. the time folder to read from
	// // 	    mesh,
	// // 	    IOobject::MUST_READ, // always gets imported, will throw an error if the field is missing
	// // 	    IOobject::NO_WRITE // will get saved automatically when the controlDict parameters will request it
	// // 	),
	// // 	mesh // initialises the field to match the size of the mesh with default (0) values
	// // );

    // volVectorField U // note that velocity is a vector field
	// (
	// 	IOobject
	// 	(
	// 	    "U",
	// 	    runTimeName_,
	// 	    mesh,
	// 	    IOobject::MUST_READ,
	// 	    IOobject::NO_WRITE
	// 	),
	// 	mesh
	// );

    // forAll(cellSets_, i)
    // {
    //     const label cell = cellSets_[i];
    //     Info << "cell ID is: " << cell << ",  "
    //          << "velocity is " << U[cell] << endl;
    // }

    

    // Info << "Cells in cellset " << cellSetName_ << ":" << endl;
    // forAll(cellSets_, i)
    // {
    //     const label cell = cellSets_[i];
    //     Info << "cell ID is: " << cell << endl;
    //     Info << mesh.C()[cell] << endl;
    // }

    // word cellZoneName_ (customDict.lookup("cellZoneName"));
    // label zoneID = mesh.cellZones().findZoneID(cellZoneName_);

    // // if (zoneID == -1)
    // // {
    // //     FatalErrorIn("yourFunctionName")
    // //         << "Cannot find cellZone " << cellSetName << endl
    // //         << "Valid cellZones are " << mesh.cellZones().names()
    // //         << exit(FatalError);
    // // }

    // const labelList& cellZones_ = mesh.cellZones()[zoneID];

    // Info << "Cells in cellzone " << cellZoneName_ << ":" << endl;
    // forAll(cellZones_, i)
    // {
    //     const label cell = cellZones_[i];
    //     Info << "cell ID is: " << cell << endl;
    //     Info << mesh.C()[cell] << endl;
    // }


    return 0;


}


// ************************************************************************* //
