

    // add snapshots in specific time interval
    if (runTime.value() >= snapshotsTime)
    {
        Info<< "runTime.value(): " << runTime.value() << endl
            << "snapshotsTime: " << snapshotsTime << endl;
        forAll(mesh.C(), cellI)
        {
            snapshotsM[cellI][snapshotsNo] = T[cellI];
            snapshotsLapM[cellI][snapshotsNo] = laplacianT[cellI];
        }

        snapshotsTime += timeInterval;
        ++ snapshotsNo;
    }    
