
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

    RectangularMatrix<scalar> nonLinearCoeffMatrix(modesNum, modesNum);
    List<scalar> tempList(mesh.C());
    dataFile = mesh.time().path()/"SVD"/"nonLinearCoeffMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < nonLinearCoeffMatrix.m(); ++row)
    {
        for (label column = 0; column < nonLinearCoeffMatrix.n(); ++column)
        {
            
            forAll(tempList, i)
            {
                tempList[i] = projModes(i, row) * 
            }

            RectangularMatrix<scalar> (projModes.subColumn(row));
            RectangularMatrix<scalar> (projModes.subColumn(row))
            
            outputFilePtr().width(16);
            outputFilePtr() << diffuTermCoeffMatrix[row][column];
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
    