
    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> diffuTermCoeffMatrix(modesNum, modesNum);
    RectangularMatrix<scalar> projModes(fieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
    RectangularMatrix<scalar> projLapModes(laplacianModesMatrix.subMatrix(0, 0, -1, modesNum));

    diffuTermCoeffMatrix = projModes.T() * projLapModes;

    // scalar diffuTermCoeffij;
    dataFile = mesh.time().path()/"SVD"/"diffuTermCoeffMatrix";
    writeMatrix(diffuTermCoeffMatrix, dataFile);

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
    writeMatrix(nonLinearCoeffM1, dataFile);

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
    writeMatrix(nonLinearCoeffM2, dataFile);
    