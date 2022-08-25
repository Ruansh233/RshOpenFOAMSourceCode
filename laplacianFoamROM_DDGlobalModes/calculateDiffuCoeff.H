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
                
        // // write subDdiffuCoeffMatrix to file
        // dataFile = mesh.time().path()/"SVD"/"subDdiffuCoeffMatrix" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDdiffuCoeffMatrix.m(); ++row)
        // {
        //     for (label column = 0; column < subDdiffuCoeffMatrix.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDdiffuCoeffMatrix[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        
        // // write subDTempCoeffMatrix to file
        // dataFile = mesh.time().path()/"SVD"/"subDTempCoeffMatrix" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDTempCoeffMatrix.m(); ++row)
        // {
        //     for (label column = 0; column < subDTempCoeffMatrix.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDTempCoeffMatrix[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        subDCoeffMatrix = SVDinv(subDTempCoeffMatrix) * subDdiffuCoeffMatrix;
        
        // write subDCoeffMatrix to file
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
    