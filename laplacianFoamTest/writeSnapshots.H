    fileName dataFile (mesh.time().path()/"SVD"/"SnapshotMatrix");
    autoPtr<OFstream> outputFilePtr;
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < snapshotsM.m(); ++row)
    {
        for (label column = 0; column < snapshotsM.n(); ++column)
        {
            outputFilePtr().width(12);
            outputFilePtr() << snapshotsM[row][column];
        }
        outputFilePtr() << endl;
    }
