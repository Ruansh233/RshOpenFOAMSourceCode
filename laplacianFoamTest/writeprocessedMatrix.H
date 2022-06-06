    autoPtr<OFstream> outputFilePtr;
    
    fileName dataFile (mesh.time().path()/"SVD"/"SnapshotMatrix");    
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
    
    dataFile = mesh.time().path()/"SVD"/"modeMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < fieldValue.U().m(); ++row)
    {
        for (label column = 0; column < fieldValue.U().n(); ++column)
        {
            outputFilePtr().width(12);
            outputFilePtr() << fieldValue.U()[row][column];
        }
        outputFilePtr() << endl;
    }

    dataFile = mesh.time().path()/"SVD"/"coeffMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < coeff.m(); ++row)
    {
        for (label column = 0; column < coeff.n(); ++column)
        {
            outputFilePtr().width(12);
            outputFilePtr() << coeff[row][column];
        }
        outputFilePtr() << endl;
    }

    dataFile = mesh.time().path()/"SVD"/"eigenVMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < fieldValue.S().size(); ++row)
    {
        outputFilePtr() << fieldValue.S()[row];
        outputFilePtr() << endl;
    }
