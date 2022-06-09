

    // add snapshots in specific time interval
    if (runTime.value() >= snapshotsTime)
    {
        // forAll(mesh.C(), cellI)
        // {
        //     snapshotsM[cellI][snapshotsNo] = T[cellI];
        // }
        

        for(label row = 0; row < snapshotsRows; ++row)
        {
            if(row < mesh.C().size())
            {
                snapshotsM[row][snapshotsNo] = T[row];
            }
            else
            {
                forAll(T.boundaryField(), patchI)
                {
                    forAll(T.boundaryField()[patchI], faceI)
                    {
                        snapshotsM[row][snapshotsNo] = T.boundaryField()[patchI][faceI];
                        ++row;
                    }
                }     
            }
        }

        snapshotsTime += timeInterval;
        ++ snapshotsNo;
    }    
