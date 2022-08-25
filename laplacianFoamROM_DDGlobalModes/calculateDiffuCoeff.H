    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> projModes(fieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
    RectangularMatrix<scalar> projLapModes(laplacianModesMatrix.subMatrix(0, 0, -1, modesNum));

    forAll(mesh.cellZones(), zoneI)
    {
        
        RectangularMatrix<scalar> subDProjModes(mesh.cellZones()[zoneI].size(), modesNum);
        RectangularMatrix<scalar> subDProjLapModes(mesh.cellZones()[zoneI].size(), modesNum);
        SquareMatrix<scalar> subDdiffuCoeffMatrix(modesNum, modesNum);
        RectangularMatrix<scalar> subDTempCoeffMatrix(modesNum, modesNum);

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
                
        // write subDdiffuCoeffMatrix to file
        dataFile = mesh.time().path()/"SVD"/"subDdiffuCoeffMatrix" + "domain" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDdiffuCoeffMatrix.m(); ++row)
        {
            for (label column = 0; column < subDdiffuCoeffMatrix.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDdiffuCoeffMatrix[row][column];
            }
            outputFilePtr() << endl;
        }

        // write subDTempCoeffMatrix to file
        dataFile = mesh.time().path()/"SVD"/"subDTempCoeffMatrix" + "domain" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDTempCoeffMatrix.m(); ++row)
        {
            for (label column = 0; column < subDTempCoeffMatrix.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDTempCoeffMatrix[row][column];
            }
            outputFilePtr() << endl;
        }

        // // write subDProjModes to file
        // dataFile = mesh.time().path()/"SVD"/"subDProjModes" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDProjModes.m(); ++row)
        // {
        //     for (label column = 0; column < subDProjModes.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDProjModes[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }
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
    