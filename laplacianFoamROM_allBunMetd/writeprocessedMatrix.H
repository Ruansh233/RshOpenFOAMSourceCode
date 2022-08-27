    Info << "write svd data\n" << endl;
    
    // check whether SVD folder exist or not
    if(!isDir(dataFile = mesh.time().path()/"SVD"))
        mkDir(dataFile);    
        
    // write SnapshotMatrix
    dataFile = mesh.time().path()/"SVD"/"SnapshotMatrix";    
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
    Info << "write Snapshot Matrix: " << snapshotsM.sizes() << endl;

    // write spatial modeMatrix
    dataFile = mesh.time().path()/"SVD"/"modeMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < fieldValue.U().m(); ++row)
    {
        for (label column = 0; column < fieldValue.U().n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << fieldValue.U()[row][column];
        }
        outputFilePtr() << endl;
    }
    Info << "write spatial mode Matrix: " << fieldValue.U().sizes() << endl;

    // RectangularMatrix<scalar> testM(fieldValue.U().n(), fieldValue.U().m());
    // dataFile = mesh.time().path()/"SVD"/"testM";
    // outputFilePtr.reset(new OFstream(dataFile));
    // testM = fieldValue.U().T();
    // // testM = multiply(fieldValue.U(), fieldValue.U().T());
    // // Info << testM;
    // for (label row = 0; row < testM.m(); ++row)
    // {
    //     for (label column = 0; column < testM.n(); ++column)
    //     {
    //         outputFilePtr().width(16);
    //         outputFilePtr() << testM[row][column];
    //     }
    //     outputFilePtr() << endl;
    // }
    // // Info << "test: \n" << fieldValue.V() * fieldValue.V().T() << endl;

    // write temporal coefficient Matrix
    RectangularMatrix<scalar> coeff(snapshotsM.n(), snapshotsM.n());
    coeff = snapshotsM.T() * fieldValue.U();

    dataFile = mesh.time().path()/"SVD"/"coeffMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < coeff.m(); ++row)
    {
        // for (label column = 0; column < coeff.n(); ++column)
        for (label column = 0; column < coeff.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << coeff[row][column];
        }
        outputFilePtr() << endl;
    }
    Info << "write temporal coefficient Matrix: " << coeff.sizes() << endl;

    // RectangularMatrix<scalar> testT(snapshotsM.m(), snapshotsM.n());
    // RectangularMatrix<scalar> testphi(fieldValue.U().subMatrix(0, 0, -1, 6));
    // RectangularMatrix<scalar> testcoeff(coeff.subMatrix(0, 0, -1, 6));
    // testT = testphi * testcoeff.T();

    // dataFile = mesh.time().path()/"SVD"/"testTMatrix";
    // outputFilePtr.reset(new OFstream(dataFile));
    // for (label row = 0; row < testT.m(); ++row)
    // {
    //     // for (label column = 0; column < coeff.n(); ++column)
    //     for (label column = 0; column < testT.n(); ++column)
    //     {
    //         outputFilePtr().width(16);
    //         outputFilePtr() << testT[row][column];
    //     }
    //     outputFilePtr() << endl;
    // }

    // Info << "testphi: " << testphi.sizes() << endl
    //      << "testcoeff: " << testcoeff.sizes() << endl;

    // write eigenvalue of covariance Matrix
    // eigenvalue is squre of the diagnoal Matrix in SVD
    dataFile = mesh.time().path()/"SVD"/"eigenVMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    // scalar sumofEigenV (0.0);
    // for (label row = 0; row < fieldValue.S().size(); ++row)
    // {
    //     sumofEigenV += pow(fieldValue.S()[row], 2);
    // }

    // for (label row = 0; row < fieldValue.S().size(); ++row)
    // {
    //     outputFilePtr() << pow(fieldValue.S()[row], 2)/sumofEigenV;
    //     outputFilePtr() << endl;
    // }
    for (label row = 0; row < fieldValue.S().size(); ++row)
    {
        outputFilePtr() << fieldValue.S()[row];
        outputFilePtr() << endl;
    }
    Info << "write eigenvalue: " << fieldValue.S().size() << endl;

    // write laplacian of modes Matrix
    dataFile = mesh.time().path()/"SVD"/"laplacianModesMatrix";
    outputFilePtr.reset(new OFstream(dataFile));

    for (label row = 0; row < laplacianModesMatrix.m(); ++row)
    {
        for (label column = 0; column < laplacianModesMatrix.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << laplacianModesMatrix[row][column];
        }
        outputFilePtr() << endl;
    }
    Info << "laplacian of modes Matrix: " << laplacianModesMatrix.sizes() << endl;

