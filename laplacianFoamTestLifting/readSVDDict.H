
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

    scalar snapshotsTime (0);
    scalar timeInterval  (svdDict.lookupOrDefault<scalar>("timeInterval", 0.1));
    label defaultSnapshotsNum (label(runTime.endTime().value()/timeInterval + 1));
    // snapshotsNum: snapshots number from dict
    label snapshotsNum (svdDict.lookupOrDefault<label>("snapshots", defaultSnapshotsNum));
    // snapshotsNo: snapshots extract in running, using to swith to new column
    label snapshotsNo (0);
    // how many mode needed
    label modesNum (svdDict.lookupOrDefault<label>("modeExtract", 5));

    
    // total number of boundary faces
    label boundaryFacesNum(0);
    forAll(mesh.boundary(), patchI)
    {
        boundaryFacesNum += mesh.boundary()[patchI].size();
    }
    // total number of snapshots rows
    // label snapshotsRows (boundaryFacesNum + mesh.C().size());
    label snapshotsRows (mesh.C().size());
    RectangularMatrix<scalar> snapshotsM(snapshotsRows, snapshotsNum);

    Info << "snapshots extract\n" << "time Interval of snapshots is: " << timeInterval << endl << endl;

    // for(label row = 0; row < snapshotsRows; ++row)
    // {
    //     if(row < mesh.C().size())
    //     {
    //         snapshotsM[row][snapshotsNo] = T[row];
    //     }
    //     else
    //     {
    //         forAll(T.boundaryField(), patchI)
    //         {
    //             forAll(T.boundaryField()[patchI], faceI)
    //             {
    //                 snapshotsM[row][snapshotsNo] = T.boundaryField()[patchI][faceI];
    //                 ++row;
    //             }
    //         }     
    //     }
    // }

    forAll(mesh.C(), cellI)
    {
        snapshotsM[cellI][snapshotsNo] = T[cellI] - liftingFunction[cellI];
    }

    snapshotsTime += timeInterval;
    ++ snapshotsNo;
