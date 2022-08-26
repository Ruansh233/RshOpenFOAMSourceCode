    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> projModes(fieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
    RectangularMatrix<scalar> projLapModes(laplacianModesMatrix.subMatrix(0, 0, -1, modesNum));

    forAll(mesh.cellZones(), zoneI)
    {
        
        RectangularMatrix<scalar> subDProjModes(mesh.cellZones()[zoneI].size(), modesNum);
        RectangularMatrix<scalar> subDProjLapModes(mesh.cellZones()[zoneI].size(), modesNum);

        RectangularMatrix<scalar> subDdiffuCoeffMatrix(modesNum, modesNum);
        RectangularMatrix<scalar> subDTempCoeffMatrix(modesNum, modesNum);
        RectangularMatrix<scalar> subDCoeffMatrix(modesNum, modesNum);

        for(label col=0; col < modesNum; ++col)
        {  
            forAll(mesh.cellZones()[zoneI], cellI)
            {
                label cell = mesh.cellZones()[zoneI][cellI];
                subDProjModes(cellI, col) = projModes(cell, col);
                subDProjLapModes(cellI, col) = projLapModes(cell, col);               
            }
        }

        subDTempCoeffMatrix = subDProjModes.T() * subDProjModes;
        subDdiffuCoeffMatrix = subDProjModes.T() * subDProjLapModes;
        subDCoeffMatrix = SVDinv(subDTempCoeffMatrix) * subDdiffuCoeffMatrix;
        
        // write subdomain ROM CoeffMatrix to file
        dataFile = mesh.time().path()/"SVD"/"subDCoeffMatrix" + "domain" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDCoeffMatrix.m(); ++row)
        {
            for (label column = 0; column < subDCoeffMatrix.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDCoeffMatrix[row][column];
            }
            outputFilePtr() << endl;
        }

        // write subdomain temporal coefficients
        RectangularMatrix<scalar> subDfieldValue(mesh.cellZones()[zoneI].size(), snapshotsM.n());
        RectangularMatrix<scalar> subDtemporalCoeff(mesh.cellZones()[zoneI].size(), modesNum);

        for(label col=0; col < fieldValueSVD.U().n(); ++col)
        {  
            forAll(mesh.cellZones()[zoneI], cellI)
            {
                label cell = mesh.cellZones()[zoneI][cellI];
                subDfieldValue(cellI, col) = snapshotsM(cell, col); 
            }
        }

        subDtemporalCoeff = subDfieldValue.T() * subDProjModes;

        // write subDtemporalCoeff to file
        dataFile = mesh.time().path()/"SVD"/"subDtemporalCoeff" + "domain" + name(zoneI);
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
        
    }

    // diffuTermCoeffMatrix = projModes.T() * projLapModes;

    // // scalar diffuTermCoeffij;
    // dataFile = mesh.time().path()/"SVD"/"diffuTermCoeffMatrix";
    // outputFilePtr.reset(new OFstream(dataFile));
    // for (label row = 0; row < diffuTermCoeffMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < diffuTermCoeffMatrix.n(); ++column)
    //     {
    //         outputFilePtr().width(16);
    //         outputFilePtr() << diffuTermCoeffMatrix[row][column];
    //     }
    //     outputFilePtr() << endl;
    // }

    // calculate inverse of a matrix
    // SquareMatrix<scalar> tt(2, 2);
    //     tt(0, 0) = 1;
    //     tt(0, 1) = 2;
    //     tt(1, 0) = 3;
    //     tt(1, 1) = 4;

    // Info<< SVDinv(tt) << endl;
    