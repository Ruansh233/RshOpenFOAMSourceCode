        // add snapshots in specific time interval
        if (runTime.value() >= snapshotsTime)
        {
            forAll(mesh.C(), cellI)
            {
                snapshotsM[cellI][snapshotsNo] = T[cellI];
            }
            snapshotsTime += timeInterval;
            ++ snapshotsNo;
        }    
