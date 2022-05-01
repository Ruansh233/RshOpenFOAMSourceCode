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
#include "wordRe.H"

int main(int argc, char *argv[])
{
    // add timeSelector
    timeSelector::addOptions();

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

    // dict Selector
    const word dictName("cellZoneDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

    // create dictionary
    IOdictionary customDict(dictIO);

    Info << customDict;

    // // Get access to a custom dictionary
    // // dictionary customDict;
    // const word dictName("cellZoneDict");

    // // Create and input-output object - this holds the path to the dict and its name
    // IOdictionary customDict
    // (
    //     IOobject
    //     (
    //         dictName, // name of the file
    //         mesh.time().system(), // path to where the file is
    //         mesh, // reference to the mesh needed by the constructor
    //         IOobject::MUST_READ // indicate that reading this dictionary is compulsory
    //     )
    // );

    // Info << "mesh.time(): " << mesh.time().name() << endl;
    // Info << "runTime.timeName(): " << runTime.timeName() << endl;

    if (mesh.cellZones().size() == 0)
    {
        FatalErrorIn("cellZoneFieldValue")
            << "There is no cellZone in this mesh"
            << exit(FatalError);

    }

    // // read variable from list
    // word runTimeName_ (customDict.lookup("readTime"));
    // word fieldName_ (customDict.lookup("fieldName"));
    List<word> cellZonesName_ (customDict.lookup("cellZoneName"));
    Info << "cellZonesName: " << cellZonesName_ << endl;
    List<word> fieldsName_ (customDict.lookup("fieldName"));

    // timeSelector 
    instantList timeDirs = timeSelector::select0(runTime, args);
    forAll(timeDirs, timei)
    {
        word runTimeName_ (timeDirs[timei].name());
        Info << "run time processed: " << runTimeName_ << endl;

        // read different field values depend on fieldName 
        // List<scalar> field_;

        fileName outputDir = mesh.time().path()/"postProcessing";
        mkDir(outputDir);
        // word fileName_(fieldsName_[0]);

        #include "createField.H"
    }



    // volVectorField U // note that velocity is a vector field
    // (
    //     IOobject
    //     (
    //         "U",
    //         runTimeName_,
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh
    // );

    // volScalarField field_(U.component(vector::Z));
    
    // // create list for writing file
    // List<scalarList> outputData_(cellZonesName_.size());

    // fileName outputDir = mesh.time().path()/"postProcessing";
    // mkDir(outputDir);
    // OFstream outputFile_(outputDir/fieldName_);
    
    // if (cellZonesName_[0] == "allZone")
    // {
    //     outputFile_.width(10);
    //     outputFile_ << "cell ID";      
    //     outputData_.resize(mesh.cellZones().size());
    //     forAll(mesh.cellZones(), zoneI)
    //     {
    //         Info << endl << "cellZone name: " << mesh.cellZones()[zoneI].name() << endl;
    //         outputFile_.width(16);
    //         outputFile_ << mesh.cellZones()[zoneI].name();

    //         forAll(mesh.cellZones()[zoneI], cellI)
    //         {
    //             const label cell = mesh.cellZones()[zoneI][cellI];
    //             outputData_[zoneI].append(field_[cell]);
    //         }            
    //     }
    //     outputFile_ << endl;
    // }
    // else
    // {
    //     outputFile_.width(10);
    //     outputFile_ << "cell ID";
    //     forAll(cellZonesName_, zoneI)
    //     {
            
    //         Info << endl << "cellZone name: " << cellZonesName_[zoneI] << endl;
    //         outputFile_.width(16);
    //         outputFile_ << cellZonesName_[zoneI];
            
    //         label zoneID = mesh.cellZones().findZoneID(cellZonesName_[zoneI]);

    //         if (zoneID == -1)
    //         {
    //             FatalErrorIn("yourFunctionName")
    //                 << "Cannot find cellZone " << cellZonesName_[zoneI] << endl
    //                 << "Valid cellZones are " << mesh.cellZones().names()
    //                 << exit(FatalError);
    //         }

    //         forAll(mesh.cellZones()[zoneID], cellI)
    //         {
    //             const label cell = mesh.cellZones()[zoneID][cellI];
    //             outputData_[zoneI].append(field_[cell]);
    //         }
    //     }
    //     outputFile_ << endl;
    // }
    

    

    // forAll(outputData_[0], dataID)
    // {
    //     word rowName_;
    //     rowName_ = "cell " + name(dataID);
    //     outputFile_.width(10);
    //     outputFile_ << rowName_;
        
    //     forAll(outputData_, blockID)
    //     {
    //         outputFile_.width(16);
    //         outputFile_ << outputData_[blockID][dataID];
    //     }
    //     outputFile_ << endl;
    // }

    Info << "\nEnd\n" << endl;

    return 0;

}


// ************************************************************************* //
