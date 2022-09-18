
    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> diffuTermCoeffMatrix(modesNum, modesNum);
    RectangularMatrix<scalar> projModes(fieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
    RectangularMatrix<scalar> projLapModes(laplacianModesMatrix.subMatrix(0, 0, -1, modesNum));

    diffuTermCoeffMatrix = projModes.T() * projLapModes;

    // scalar diffuTermCoeffij;
    dataFile = mesh.time().path()/"SVD"/"diffuTermCoeffMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < diffuTermCoeffMatrix.m(); ++row)
    {
        for (label column = 0; column < diffuTermCoeffMatrix.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << diffuTermCoeffMatrix[row][column];
        }
        outputFilePtr() << endl;
    }

    scalar tensorM = modesNum * modesNum;
    RectangularMatrix<scalar> nonLinearCoeffM(tensorM, modesNum);
    List<scalar> tempList(mesh.C().size());
    for (label i = 0; i < nonLinearCoeffM.n(); ++i)
    {
        for (label j = 0; j < nonLinearCoeffM.n(); ++j)
        {
            for (label k = 0; k < nonLinearCoeffM.n(); ++k)
            {            
                forAll(tempList, listI)
                {
                    tempList[listI] = fieldValueSVD.U()(listI, i) * fieldValueSVD.U()(listI, j) * laplacianModesMatrix(listI, k);
                }
                nonLinearCoeffM(i * modesNum + j, k) = gSum(tempList);
                Info << "j, k: " << i * modesNum + j << ", " << k << ": " 
                     << nonLinearCoeffM(i * modesNum + j, k) << endl;
            }
        }
    }

    dataFile = mesh.time().path()/"SVD"/"nonLinearCoeffM";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < nonLinearCoeffM.m(); ++row)
    {
        for (label column = 0; column < nonLinearCoeffM.n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << nonLinearCoeffM[row][column];
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
    