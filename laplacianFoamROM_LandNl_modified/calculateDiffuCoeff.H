
    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> diffuTermCoeffMatrix(modesNum, modesNum);
    RectangularMatrix<scalar> projModes(fieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
    RectangularMatrix<scalar> projLapModes(laplacianModesMatrix.subMatrix(0, 0, -1, modesNum));

    diffuTermCoeffMatrix = projModes.T() * projLapModes;

    // scalar diffuTermCoeffij;
    dataFile = mesh.time().path()/"SVD"/"diffuTermCoeffMatrix";
    writeMatrix(diffuTermCoeffMatrix, dataFile);
    
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

    scalar tensorM = modesNum * modesNum;
    RectangularMatrix<scalar> nonLinearCoeffM1(tensorM, modesNum);
    List<scalar> tempList(mesh.C().size());
    for (label i = 0; i < nonLinearCoeffM1.n(); ++i)
    {
        for (label j = 0; j < nonLinearCoeffM1.n(); ++j)
        {
            for (label k = 0; k < nonLinearCoeffM1.n(); ++k)
            {            
                forAll(tempList, listI)
                {
                    tempList[listI] = fieldValueSVD.U()(listI, i) * fieldValueSVD.U()(listI, j) * laplacianModesMatrix(listI, k);
                }
                nonLinearCoeffM1(i * modesNum + j, k) = gSum(tempList);
            }
        }
    }

    dataFile = mesh.time().path()/"SVD"/"nonLinearCoeffM1";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < nonLinearCoeffM1.m(); ++row)
    {
        for (label column = 0; column < nonLinearCoeffM1.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << nonLinearCoeffM1[row][column];
        }
        outputFilePtr() << endl;
    }

    RectangularMatrix<scalar> nonLinearCoeffM2(tensorM, modesNum);
    for (label i = 0; i < nonLinearCoeffM2.n(); ++i)
    {
        for (label j = 0; j < nonLinearCoeffM2.n(); ++j)
        {
            for (label k = 0; k < nonLinearCoeffM2.n(); ++k)
            {            
                forAll(tempList, listI)
                {
                    tempList[listI] = fieldValueSVD.U()(listI, i) * (gradModesMatrix(listI, j) & gradModesMatrix(listI, k));
                }
                nonLinearCoeffM2(i * modesNum + j, k) = gSum(tempList);
            }
        }
    }

    dataFile = mesh.time().path()/"SVD"/"nonLinearCoeffM2";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < nonLinearCoeffM2.m(); ++row)
    {
        for (label column = 0; column < nonLinearCoeffM2.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << nonLinearCoeffM2[row][column];
        }
        outputFilePtr() << endl;
    }

    //     outputFilePtr() << endl;
    // }

    // scalar domainVolume (gSum(mesh.V()));

    // for (label row = 0; row < diffuTermCoeffMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < diffuTermCoeffMatrix.n(); ++column)
    //     {
    //         diffuTermCoeffij = 0;
    //         forAll(mesh.C(), cellI)
    //         {
    //             diffuTermCoeffij += 
    //             fieldValue.U()[cellI][row] * laplacianModesMatrix[cellI][column] * 
    //             mesh.V()[cellI];
    //         }

    //         diffuTermCoeffMatrix[row][column] = diffuTermCoeffij;    

    //         outputFilePtr().width(16);
    //         outputFilePtr() << diffuTermCoeffMatrix[row][column];
    //     }

    //     outputFilePtr() << endl;
    // }
    