    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> projModes(fieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
    RectangularMatrix<scalar> projLapModes(laplacianModesMatrix.subMatrix(0, 0, -1, modesNum));

    forAll(mesh.cellZones(), zoneI)
    {
        
        RectangularMatrix<scalar> subDProjModes(mesh.cellZones()[zoneI].size(), modesNum);
        RectangularMatrix<scalar> subDProjLapModes(mesh.cellZones()[zoneI].size(), modesNum);

        // coefficients of diffusion term in laplacian function
        RectangularMatrix<scalar> subDdiffuCoeffMatrix(modesNum, modesNum);
        // coefficients of transient term in laplacian function
        RectangularMatrix<scalar> subDTransCoeffMatrix(modesNum, modesNum);
        // coefficients of matrix form
        RectangularMatrix<scalar> subDROMCoeffMatrix(modesNum, modesNum);

        // assign value to proj matrix of each subdomain
        for(label col=0; col < modesNum; ++col)
        {  
            forAll(mesh.cellZones()[zoneI], cellI)
            {
                label cell = mesh.cellZones()[zoneI][cellI];
                subDProjModes(cellI, col) = projModes(cell, col);
                subDProjLapModes(cellI, col) = projLapModes(cell, col);               
            }
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

        // // write subDProjLapModes to file
        // dataFile = mesh.time().path()/"SVD"/"subDProjLapModes" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDProjLapModes.m(); ++row)
        // {
        //     for (label column = 0; column < subDProjLapModes.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDProjLapModes[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // phi.j * phi.i, transient term
        subDTransCoeffMatrix = subDProjModes.T() * subDProjModes;
        // laplacian(phi.j) * phi.i， diffusion term
        subDdiffuCoeffMatrix = subDProjModes.T() * subDProjLapModes;

        // write subDTransCoeffMatrix to file
        dataFile = mesh.time().path()/"SVD"/"subDTransCoeffMatrix" + "domain" + name(zoneI);
        outputFilePtr.reset(new OFstream(dataFile));
        for (label row = 0; row < subDTransCoeffMatrix.m(); ++row)
        {
            for (label column = 0; column < subDTransCoeffMatrix.n(); ++column)
            {
                outputFilePtr().width(16);
                outputFilePtr() << subDTransCoeffMatrix[row][column];
            }
            outputFilePtr() << endl;
        }

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

        // subDROMCoeffMatrix = SVDinv(subDTransCoeffMatrix) * subDdiffuCoeffMatrix;
        
        // // write subdomain ROM CoeffMatrix to file
        // dataFile = mesh.time().path()/"SVD"/"subDROMCoeffMatrix" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDROMCoeffMatrix.m(); ++row)
        // {
        //     for (label column = 0; column < subDROMCoeffMatrix.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDROMCoeffMatrix[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // // write SVDinv(subDTransCoeffMatrix) to file
        // dataFile = mesh.time().path()/"SVD"/"invsubDROMCoeffMatrix" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < SVDinv(subDTransCoeffMatrix).m(); ++row)
        // {
        //     for (label column = 0; column < SVDinv(subDTransCoeffMatrix).n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << SVDinv(subDTransCoeffMatrix)[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // // test the SVDinv function
        // RectangularMatrix<scalar> subDdiffuCoeffMatrixTest(modesNum, modesNum);
        // subDdiffuCoeffMatrixTest = subDTransCoeffMatrix * subDROMCoeffMatrix;

        // // write subDdiffuCoeffMatrixTest to file
        // dataFile = mesh.time().path()/"SVD"/"subDdiffuCoeffMatrixTest" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDdiffuCoeffMatrixTest.m(); ++row)
        // {
        //     for (label column = 0; column < subDdiffuCoeffMatrixTest.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDdiffuCoeffMatrixTest[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

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

        // // write subdomain temporal coefficients
        // // subDfieldValue, subdomain field value
        // // subDtemporalCoeff, subdomain temporal coefficients
        // RectangularMatrix<scalar> subDfieldValue(mesh.cellZones()[zoneI].size(), snapshotsM.n());
        // for(label col=0; col < fieldValueSVD.U().n(); ++col)
        // {  
        //     forAll(mesh.cellZones()[zoneI], cellI)
        //     {
        //         label cell = mesh.cellZones()[zoneI][cellI];
        //         subDfieldValue(cellI, col) = snapshotsM(cell, col); 
        //     }
        // }

        // // write subDfieldValue to file
        // dataFile = mesh.time().path()/"SVD"/"subDfieldValue" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDfieldValue.m(); ++row)
        // {
        //     for (label column = 0; column < subDfieldValue.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDfieldValue[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // // A=T.T*phi*[(phi.T*phi)].T
        // RectangularMatrix<scalar> subDtemporalCoeff(mesh.cellZones()[zoneI].size(), modesNum);
        // subDtemporalCoeff = subDfieldValue.T() * subDProjModes * 
        //                     SVDinv(subDProjModes.T() * subDProjModes).T();

        // // write subDtemporalCoeff to file
        // dataFile = mesh.time().path()/"SVD"/"subDtemporalCoeff" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDtemporalCoeff.m(); ++row)
        // {
        //     for (label column = 0; column < subDtemporalCoeff.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDtemporalCoeff[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // // temporal coefficient calculation test
        // RectangularMatrix<scalar> coeffTest(snapshotsM.n(), snapshotsM.n());    
        // for (label column = 0; column < coeffTest.n(); ++column)
        // {
        //     RectangularMatrix<scalar> coeffTestTemp(fieldValueSVD.V().subColumn(column, 0, -1));
        //     for (label row = 0; row < coeffTest.m(); ++row)
        //     {
        //         coeffTest(row, column) = fieldValueSVD.S()[column] * coeffTestTemp(row, 0);
        //     }
        // }

        // dataFile = mesh.time().path()/"SVD"/"coeffMatrixTest";
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < coeffTest.m(); ++row)
        // {
        //     for (label column = 0; column < coeffTest.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << coeffTest[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // RectangularMatrix<scalar> subDfieldValueTest2(mesh.cellZones()[zoneI].size(), snapshotsM.n());
        // subDfieldValueTest2 = subDProjModes * coeffTest.T();

        // // write test subDfieldValue to file
        // dataFile = mesh.time().path()/"SVD"/"subDfieldValueTest2" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDfieldValueTest2.m(); ++row)
        // {
        //     for (label column = 0; column < subDfieldValueTest2.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDfieldValueTest2[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }  



        // RectangularMatrix<scalar> subDfieldValueTest(mesh.cellZones()[zoneI].size(), snapshotsM.n());
        // subDfieldValueTest = subDProjModes * subDtemporalCoeff.T();

        // // write test subDfieldValue to file
        // dataFile = mesh.time().path()/"SVD"/"subDfieldValueTest" + "domain" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDfieldValueTest.m(); ++row)
        // {
        //     for (label column = 0; column < subDfieldValueTest.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDfieldValueTest[row][column];
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

    