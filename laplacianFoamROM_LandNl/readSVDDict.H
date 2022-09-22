
    const word dictName("svdDict");

    IOdictionary svdDict
    (
        IOobject
        (
            dictName, // name of the file
            mesh.time().system(), // path to where the file is
            mesh, // reference to the mesh needed by the constructor
            IOobject::MUST_READ, // indicate that reading this dictionary is compulsory
            IOobject::NO_WRITE
        )

    );

    bool ROMorNot (svdDict.lookupOrDefault<bool>("ROMorNot", false));
    
    scalar snapshotsTime (0);
    scalar timeInterval  (svdDict.lookupOrDefault<scalar>("timeInterval", 0.1));
    label snapshotsNum (label(runTime.endTime().value()/timeInterval + 1));
    // snapshotsNo: snapshots extract in running, using to swith to new column
    label snapshotsNo (0);
    // how many mode needed
    label modesNum (svdDict.lookupOrDefault<label>("modeExtract", 5));
    label snapshotsRows (mesh.C().size());
    RectangularMatrix<scalar> snapshotsM(snapshotsRows, snapshotsNum);
    RectangularMatrix<scalar> snapshotsLapM(snapshotsRows, snapshotsNum);
    
    if(ROMorNot)
    {       
        // // total number of boundary faces
        // label boundaryFacesNum(0);
        // forAll(mesh.boundary(), patchI)
        // {
        //     boundaryFacesNum += mesh.boundary()[patchI].size();
        // }
        // // total number of snapshots rows
        // // label snapshotsRows (boundaryFacesNum + mesh.C().size());

        Info << "snapshots extract\n" << "time Interval of snapshots is: " << timeInterval << endl << endl;

        forAll(mesh.C(), cellI)
        {
            snapshotsM[cellI][snapshotsNo] = T[cellI];
        }

        snapshotsTime += timeInterval;
        ++ snapshotsNo;
    }

    
