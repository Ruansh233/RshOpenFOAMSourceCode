    Info << "write svd data\n" << endl;
        
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
    for (label row = 0; row < fieldValueSVD.U().m(); ++row)
    {
        for (label column = 0; column < fieldValueSVD.U().n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << fieldValueSVD.U()[row][column];
        }
        outputFilePtr() << endl;
    }
    Info << "write spatial mode Matrix: " << fieldValueSVD.U().sizes() << endl;

    RectangularMatrix<scalar> testM(fieldValueSVD.U().n(), fieldValueSVD.U().m());
    dataFile = mesh.time().path()/"SVD"/"testM";
    outputFilePtr.reset(new OFstream(dataFile));
    testM = fieldValueSVD.V();
    // testM = multiply(fieldValueSVD.U(), fieldValueSVD.U().T());
    // Info << testM;
    for (label row = 0; row < testM.m(); ++row)
    {
        for (label column = 0; column < testM.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << testM[row][column];
        }
        outputFilePtr() << endl;
    }

    // write temporal coefficient Matrix
    RectangularMatrix<scalar> coeff(snapshotsM.n(), snapshotsM.n());
    coeff = snapshotsM.T() * fieldValueSVD.U();

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

    // write eigenvalue of covariance Matrix
    // eigenvalue is squre of the diagnoal Matrix in SVD
    dataFile = mesh.time().path()/"SVD"/"eigenVMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    // scalar sumofEigenV (0.0);
    // for (label row = 0; row < fieldValueSVD.S().size(); ++row)
    // {
    //     sumofEigenV += pow(fieldValueSVD.S()[row], 2);
    // }

    // for (label row = 0; row < fieldValueSVD.S().size(); ++row)
    // {
    //     outputFilePtr() << pow(fieldValueSVD.S()[row], 2)/sumofEigenV;
    //     outputFilePtr() << endl;
    // }
    for (label row = 0; row < fieldValueSVD.S().size(); ++row)
    {
        outputFilePtr() << fieldValueSVD.S()[row];
        outputFilePtr() << endl;
    }
    Info << "write eigenvalue: " << fieldValueSVD.S().size() << endl;

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

