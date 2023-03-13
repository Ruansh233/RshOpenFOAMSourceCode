    Info << "write svd data\n" << endl;
    
    // write SnapshotMatrix
    dataFile = mesh.time().path()/"SVD"/"SnapshotMatrix";    
    writeMatrix(snapshotsM, dataFile);
    Info << "write Snapshot Matrix: " << snapshotsM.sizes() << endl;

    // write spatial modeMatrix
    dataFile = mesh.time().path()/"SVD"/"modeMatrix";
    writeMatrix(fieldValueSVD.U(), dataFile);
    Info << "write spatial mode Matrix: " << fieldValueSVD.U().sizes() << endl;

    // write temporal coefficient Matrix
    RectangularMatrix<scalar> coeff(snapshotsM.n(), snapshotsM.n());
    coeff = snapshotsM.T() * fieldValueSVD.U();

    dataFile = mesh.time().path()/"SVD"/"coeffMatrix";
    writeMatrix(coeff, dataFile);
    Info << "write temporal coefficient Matrix: " << coeff.sizes() << endl;

    // write eigenvalue of covariance Matrix
    // eigenvalue is squre of the diagnoal Matrix in SVD
    dataFile = mesh.time().path()/"SVD"/"eigenVMatrix";
    writeMatrix(fieldValueSVD.S(), dataFile);
    Info << "write eigenvalue: " << fieldValueSVD.S().size() << endl;

    // write laplacian of modes Matrix
    dataFile = mesh.time().path()/"SVD"/"laplacianModesMatrix";
    writeMatrix(laplacianModesMatrix, dataFile);
    Info << "laplacian of modes Matrix: " << laplacianModesMatrix.sizes() << endl;
    