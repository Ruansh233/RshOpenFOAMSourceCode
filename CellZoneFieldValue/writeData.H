// Info << "hello";

// List<scalarList> outputData_(cellZonesName_.size());    

// // Rsh, test Matrix
// #include "RectangularMatrix.H"

// label Mrows(10);
// label Ncolumns(20);

// // RectangularMatrix<scalar> testM(Mrows,Ncolumns);


// testM[0][0] = 100.;

// Info << "testM: " << endl << testM;

List<scalarField> outputData_(cellZonesName_.size());    
OFstream outputFile_(outputDir/fileName_);

if (cellZonesName_[0] == "allZone")
{
    // outputFile_.width(10);
    // outputFile_ << "cell ID";   

    outputData_.resize(mesh.cellZones().size());
    forAll(mesh.cellZones(), zoneI)
    {
        Info << endl << "cellZone name: " << mesh.cellZones()[zoneI].name() << endl;

        outputFile_.width(16);
        outputFile_ << mesh.cellZones()[zoneI].name();
        outputData_[zoneI].resize(mesh.cellZones()[zoneI].size());

        forAll(mesh.cellZones()[zoneI], cellI)
        {
            // mesh.cellZones()[zoneI][cellI] give the global cell number in the mesh
            const label cell = mesh.cellZones()[zoneI][cellI];
            // outputData_[zoneI].append(field_[cell]);
            
            outputData_[zoneI][cellI] = field_[cell];
        }            
    }
    outputFile_ << endl;
}

else if (cellZonesName_[0].find(".*") != string::npos)
{
    // Info << "cellZonesName: " << cellZonesName_[0] << endl;
    
    wordRe cellZonesNameRe_(cellZonesName_[0]);
    cellZonesNameRe_.compile ();

    List<label> matchZoneID(mesh.cellZones().size());
    label countNumber(0);

    forAll(mesh.cellZones(), ZoneI)
    {
        
        if(cellZonesNameRe_.match(mesh.cellZones()[ZoneI].name()))
        {
            // Info << endl << "cellZone name: " << mesh.cellZones()[zoneI].name() << endl;
            
            matchZoneID[countNumber] = ZoneI;
            countNumber += 1;
        }
    }

    matchZoneID.resize(countNumber);
    outputData_.resize(matchZoneID.size());

    // Info << matchZoneID << endl;

    forAll(matchZoneID, i)
    {
        label zoneI (matchZoneID[i]);

        Info << endl << "cellZone name: " << mesh.cellZones()[zoneI].name() << endl;

        outputFile_.width(16);
        outputFile_ << mesh.cellZones()[zoneI].name();
        outputData_[i].resize(mesh.cellZones()[zoneI].size());

        // Info << mesh.cellZones()[matchZoneID[0]].size() << endl;

        forAll(mesh.cellZones()[zoneI], cellI)
        {
            // mesh.cellZones()[zoneI][cellI] give the global cell number in the mesh
            const label cell = mesh.cellZones()[zoneI][cellI];
            // Info << cell << endl;
            // outputData_[zoneI].append(field_[cell]);
     
            outputData_[i][cellI] = field_[cell];
        }            
    }

    outputFile_ << endl;
}

else
{
    // outputFile_.width(10);
    // outputFile_ << "cell ID";

    // Info << "cellZonesName_ " << cellZonesName_ << endl;

    forAll(cellZonesName_, zoneI)
    {
        
        Info << endl << "cellZone name: " << cellZonesName_[zoneI] << endl;

        outputFile_.width(16);
        outputFile_ << cellZonesName_[zoneI];
        
        label zoneID = mesh.cellZones().findZoneID(cellZonesName_[zoneI]);
        outputData_[zoneI].resize(mesh.cellZones()[zoneID].size());

        // if (zoneID == -1)
        // {
        //     FatalErrorIn("yourFunctionName")
        //         << "Cannot find cellZone " << cellZonesName_[zoneI] << endl
        //         << "Valid cellZones are " << mesh.cellZones().names()
        //         << exit(FatalError);
        // }

        forAll(mesh.cellZones()[zoneID], cellI)
        {
            const label cell = mesh.cellZones()[zoneID][cellI];
            // outputData_[zoneI].append(field_[cell]);
            
            outputData_[zoneI][cellI] = field_[cell];
        }
    }

    outputFile_ << endl;
}

forAll(outputData_[0], dataID)
{
    // word rowName_;
    // rowName_ = "cell " + name(dataID);
    // outputFile_.width(10);
    // outputFile_ << rowName_;
    
    forAll(outputData_, blockID)
    {
        outputFile_.width(16);
        outputFile_ << outputData_[blockID][dataID];
    }
    outputFile_ << endl;
}

