    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> subDFieldValue;
    RectangularMatrix<scalar> subDLapFieldValue;

    forAll(mesh.cellZones(), zoneI)
    {
        RectangularMatrix<scalar> subDFieldValue(mesh.cellZones()[zoneI].size(), snapshotsNo);
        RectangularMatrix<scalar> subDLapFieldValue(mesh.cellZones()[zoneI].size(), snapshotsNo);

        // assign value to proj matrix of each subdomain
        for(label col=0; col < snapshotsNo; ++col)
        {  
            forAll(mesh.cellZones()[zoneI], cellI)
            {
                label cell = mesh.cellZones()[zoneI][cellI];
                subDFieldValue(cellI, col) = snapshotsM(cell, col);
                subDLapFieldValue(cellI, col) = snapshotsLapM(cell, col);               
            }
        }

        SVD subDfieldValueSVD(subDFieldValue); 

        // mode matrix of subdomains
        RectangularMatrix<scalar> subDModesM(mesh.cellZones()[zoneI].size(), modesNum);
        subDModesM = RectangularMatrix<scalar> (subDfieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
        // get laplacian of modes for subdomains 
        RectangularMatrix<scalar> subDLapModesM(mesh.cellZones()[zoneI].size(), modesNum);
        subDLapModesM = subDLapFieldValue * SVDinv(subDfieldValueSVD.U().T() * subDFieldValue);

        // coefficients of matrix form
        RectangularMatrix<scalar> subDROMCoeffMatrix(modesNum, modesNum);
        subDROMCoeffMatrix = subDModesM.T() * subDLapModesM;
      
        // write modes Matrix
        dataFile = mesh.time().path()/"SVD"/"subDModesM" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));

        for (label row = 0; row < subDModesM.m(); ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDModesM[row][column];
            }
            outputFilePtr() << endl;
        }
        
        // write laplacian of modes Matrix
        dataFile = mesh.time().path()/"SVD"/"subDLapModesM" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));

        for (label row = 0; row < subDLapModesM.m(); ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDLapModesM[row][column];
            }
            outputFilePtr() << endl;
        }

        // write ROM equations coefficients matrix
        dataFile = mesh.time().path()/"SVD"/"subDROMCoeffMatrix" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));

        for (label row = 0; row < subDROMCoeffMatrix.m(); ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDROMCoeffMatrix[row][column];
            }
            outputFilePtr() << endl;
        }

        // write snapshots matrix
        dataFile = mesh.time().path()/"SVD"/"subDFieldValue" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDFieldValue.m(); ++row)
        {
            for (label column = 0; column < subDFieldValue.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDFieldValue[row][column];
            }
            outputFilePtr() << endl;
        }

        // write laplacian of snapshots matrix
        dataFile = mesh.time().path()/"SVD"/"subDLapFieldValue" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDLapFieldValue.m(); ++row)
        {
            for (label column = 0; column < subDLapFieldValue.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDLapFieldValue[row][column];
            }
            outputFilePtr() << endl;
        }

        // write subDtemporalCoeff matrix
        RectangularMatrix<scalar> subDtemporalCoeff(subDFieldValue.n(), modesNum);
        subDtemporalCoeff = subDFieldValue.T() * subDModesM;
        dataFile = mesh.time().path()/"SVD"/"temporalCoeff" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDtemporalCoeff.m(); ++row)
        {
            for (label column = 0; column < subDtemporalCoeff.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDtemporalCoeff[row][column];
            }
            outputFilePtr() << endl;
        }

        // write eigen value matrix
        dataFile = mesh.time().path()/"SVD"/"eigenVMatrix";
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDfieldValueSVD.S().size(); ++row)
        {
            outputFilePtr() << subDfieldValueSVD.S()[row];
            outputFilePtr() << endl;
        }
        Info << "write eigenvalue: " << subDfieldValueSVD.S().size() << endl;

    }
