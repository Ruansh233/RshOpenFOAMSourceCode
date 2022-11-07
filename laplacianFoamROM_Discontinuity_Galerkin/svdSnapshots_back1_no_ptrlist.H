    // calculate and write matrix element value of diffusion term
    label elementCellSize (mesh.cellZones()[0].size());
    label snapshotsNum (mesh.cellZones().size());
    List<fileName> modeNames (modesNum);
    RectangularMatrix<scalar> cellFieldValue(elementCellSize, snapshotsNum);
    RectangularMatrix<scalar> faceFieldValue;

    forAll(mesh.cellZones(), zoneI)
    {
        forAll(mesh.cellZones()[zoneI], cellI)
        {
            label cell = mesh.cellZones()[zoneI][cellI];
            cellFieldValue(cellI, zoneI) = T[cell];
        }
    }

    SVD fieldValueSVD(cellFieldValue); 
    RectangularMatrix<scalar> coefficientsM(cellFieldValue.T() * fieldValueSVD.U());
    RectangularMatrix<scalar> invCoefficientsM(SVDinv(coefficientsM.T()));

    // list of boundary type wordRe
    // List<label> boundaryWordRe ({".*block.*_in", ".*block.*_out", ".*block.*_wall.*"})
    // the list word su
    List<word> boundaryWordRe ({"_in", "_out"});
    List<label> matchPatchID(mesh.boundaryMesh().size());
    wordRe matchPatch;
    label countNumber(0);

    forAll(boundaryWordRe, reI)
    {
        forAll(mesh.cellZones(), zoneI)
        {
            // Get index of patch
            // the outlet of whole domain would be put in front of all interfaces
            // --- so name(zoneI + 1) is used
            matchPatch = ".*block" + name(zoneI + 1) + boundaryWordRe[reI];
            matchPatch.compile ();

            forAll(mesh.boundaryMesh(), patchI)
            {
                if(matchPatch.match(mesh.boundary()[patchI].name()))
                {
                    matchPatchID[countNumber] = patchI;
                    ++countNumber;

                    break;
                }
            }
        }
    }
    matchPatchID.resize(countNumber);

    Info << "matchPatch: " << matchPatchID << endl;

    PtrList<RectangularMatrix<scalar>> boundarySnapshotsMList;
    PtrList<RectangularMatrix<scalar>> boundaryModesMList;
    RectangularMatrix<scalar> patchSnapshots(mesh.boundary()[matchPatchID[0]].size(), snapshotsNum);
    label patchType(0);

    forAll(matchPatchID, patchID)
    {
        forAll(mesh.boundary()[matchPatchID[patchID]], faceI)
        {
            patchSnapshots(faceI, patchID - patchType) = T.boundaryField()[matchPatchID[patchID]][faceI];
        }

        if((patchID + 1) % 5 == 0)
        {
            patchType += 5;
            boundarySnapshotsMList.append(patchSnapshots.clone());
            boundaryModesMList.append(RectangularMatrix<scalar> (patchSnapshots * invCoefficientsM).clone());

            // resides for each patchType except last patchType
            if(patchID + 1 < matchPatchID.size())
                patchSnapshots.resize(mesh.boundary()[matchPatchID[patchID + 1]].size(), snapshotsNum);
        }
    }


    forAll(boundarySnapshotsMList, listI)
    {
        dataFile = mesh.time().path()/"SVD"/"boundarySnapshots" + name(listI);
        writeMatrix(boundarySnapshotsMList[listI], dataFile);

        dataFile = mesh.time().path()/"SVD"/"boundaryModes" + name(listI);
        writeMatrix(boundaryModesMList[listI], dataFile);
    }

    dataFile = mesh.time().path()/"SVD"/"boundaryPatchID";
    writeMatrix(matchPatchID, dataFile);


    // // create fvmesh object of a single element, 
    // // but the object still contains data of whole mesh
    // fvMesh  refElementMesh
    // (
    //     IOobject
    //     (
    //         polyMesh::defaultRegion,
    //         mesh.time().path()/"refMesh",
    //         runTime,
    //         IOobject::MUST_READ
    //     ),
    //     false
    // );

    // // create mode field by copying T
    // volScalarField fieldValueMode
    // (
    //     IOobject
    //     (
    //         "TmdoeTest",
    //         mesh.time().timeName(),
    //         refElementMesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     refElementMesh,
    //     dimless
    // );

    // forAll(mesh.cellZones()[0], cellI)
    // {   
    //     label cell (mesh.cellZones()[0][cellI]);
    //     fieldValueMode[cell] = fieldValueSVD.U()(cellI, 0);
    // }
    // fieldValueMode.write();


    // // mode field value name
    // for (label iname = 0; iname < modesNum; ++iname)
    // {
    //     modeNames[iname] = "Tmode" + name(iname);
    // }

    // // set mode field value
    // forAll(modeNames, No_)
    // {
    //     // create mode field by copying T
    //     volScalarField fieldValueMode
    //     (
    //         IOobject
    //         (
    //             modeNames[No_],
    //             mesh.time().timeName(),
    //             mesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         T
    //     );

    //     // set cell value
    //     forAll(mesh.cellZones(), zoneI)
    //     {
    //         forAll(mesh.cellZones()[zoneI], cellI)
    //         {   
    //             label cell (mesh.cellZones()[zoneI][cellI]);
    //             fieldValueMode[cell] = fieldValueSVD.U()(cellI, No_);
    //         }
    //     }

    //     // set patch value
    //     patchType = 0;
    //     forAll(matchPatchID, patchID)
    //     {
    //         forAll(fieldValueMode.boundaryField()[matchPatchID[patchID]], faceI)
    //         {
    //             fieldValueMode.boundaryFieldRef()[matchPatchID[patchID]][faceI] = boundaryModesMList[patchType](faceI, No_);
    //         }
    //         if((patchID + 1) % 5 == 0)
    //         {
    //             patchType += 1;
    //         }
    //     }

    //     fieldValueMode.write();
    // }

    // // create mode field by copying T
    // volScalarField fieldValueModeTest
    // (
    //     IOobject
    //     (
    //         modeNames[0],
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );

    // // create mode field by copying T
    // volVectorField fieldValueModegrad
    // (
    //     IOobject
    //     (
    //         "grad" + modeNames[0],
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::grad(fieldValueModeTest)
    // );
    // fieldValueModegrad.write();

    // // create mode field by copying T
    // volScalarField fieldValueModelap
    // (
    //     IOobject
    //     (
    //         "lap" + modeNames[0],
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::laplacian(fieldValueModeTest)
    // );
    // fieldValueModelap.write();

    if(!isDir(mesh.time().path()/"SVD"))
        mkDir(mesh.time().path()/"SVD");

    // write snapshots Matrix
    dataFile = mesh.time().path()/"SVD"/"snapshotsM";
    writeMatrix(cellFieldValue, dataFile);

    // write modes Matrix
    dataFile = mesh.time().path()/"SVD"/"modesM";
    writeMatrix(fieldValueSVD.U(), dataFile);

    // write eigen value Matrix
    dataFile = mesh.time().path()/"SVD"/"eigenValueM";
    writeMatrix(fieldValueSVD.S(), dataFile);

    // write eigen value Matrix
    dataFile = mesh.time().path()/"SVD"/"coefficientsM";
    writeMatrix(coefficientsM, dataFile);


        // // mode matrix of subdomains
        // RectangularMatrix<scalar> subDModesM(mesh.cellZones()[zoneI].size(), modesNum);
        // subDModesM = RectangularMatrix<scalar> (subDfieldValueSVD.U().subMatrix(0, 0, -1, modesNum));
        // // get laplacian of modes for subdomains 
        // RectangularMatrix<scalar> subDLapModesM(mesh.cellZones()[zoneI].size(), modesNum);
        // subDLapModesM = subDLapFieldValue * SVDinv(subDfieldValueSVD.U().T() * subDFieldValue);

        // // // check whether the svd is correct
        // // RectangularMatrix<scalar> subDFieldValueTest(mesh.cellZones()[zoneI].size(), snapshotsNo, 0);
        // // RectangularMatrix<scalar> subDLapFieldValueTest(mesh.cellZones()[zoneI].size(), snapshotsNo, 0);

        // // for (label column = 0; column < modesNum; ++column)
        // // {
        // //     subDFieldValueTest = subDFieldValueTest + subDfieldValueSVD.S()[column] *
        // //                         RectangularMatrix<scalar> (subDModesM.subColumn(column, 0, -1)) *
        // //                         RectangularMatrix<scalar> (subDfieldValueSVD.V().T().subRow(column, 0, -1));
        // // }

        // dataFile = mesh.time().path()/"SVD"/"subDFieldValue" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDFieldValue.m(); ++row)
        // {
        //     for (label column = 0; column < subDFieldValue.n(); ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDFieldValue(row, column);
        //     }
        //     outputFilePtr() << endl;
        // }

        // // for (label column = 0; column < modesNum; ++column)
        // // {
        // //     subDLapFieldValueTest = subDLapFieldValueTest + subDfieldValueSVD.S()[column] *
        // //                         RectangularMatrix<scalar> (subDLapModesM.subColumn(column, 0, -1)) *
        // //                         RectangularMatrix<scalar> (subDfieldValueSVD.V().T().subRow(column, 0, -1));
        // // }

        // // dataFile = mesh.time().path()/"SVD"/"subDLapFieldValueTest" + name(zoneI);
        // // outputFilePtr.reset(new OFstream(dataFile));
        // // for (label row = 0; row < subDLapFieldValueTest.m(); ++row)
        // // {
        // //     for (label column = 0; column < subDLapFieldValueTest.n(); ++column)
        // //     {
        // //         outputFilePtr().width(16);
        // //         outputFilePtr() << subDLapFieldValueTest(row, column);
        // //     }
        // //     outputFilePtr() << endl;
        // // }

        // // coefficients of matrix form
        // RectangularMatrix<scalar> subDROMCoeffMatrix(modesNum, modesNum);
        // subDROMCoeffMatrix = subDModesM.T() * subDLapModesM;
      
        // // write modes Matrix
        // dataFile = mesh.time().path()/"SVD"/"subDModesM" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));

        // for (label row = 0; row < subDModesM.m(); ++row)
        // {
        //     for (label column = 0; column < modesNum; ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDModesM[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }
        
        // // write laplacian of modes Matrix
        // dataFile = mesh.time().path()/"SVD"/"subDLapModesM" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));

        // for (label row = 0; row < subDLapModesM.m(); ++row)
        // {
        //     for (label column = 0; column < modesNum; ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDLapModesM[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // // write ROM equations coefficients matrix
        // dataFile = mesh.time().path()/"SVD"/"subDROMCoeffMatrix" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));

        // for (label row = 0; row < subDROMCoeffMatrix.m(); ++row)
        // {
        //     for (label column = 0; column < modesNum; ++column)
        //     {
        //         outputFilePtr().width(16);
        //         outputFilePtr() << subDROMCoeffMatrix[row][column];
        //     }
        //     outputFilePtr() << endl;
        // }

        // // // write snapshots matrix
        // // dataFile = mesh.time().path()/"SVD"/"subDFieldValue" + name(zoneI);
        // // outputFilePtr.reset(new OFstream(dataFile));
        // // for (label row = 0; row < subDFieldValue.m(); ++row)
        // // {
        // //     for (label column = 0; column < subDFieldValue.n(); ++column)
        // //     {
        // //         outputFilePtr().width(16);
        // //         outputFilePtr() << subDFieldValue[row][column];
        // //     }
        // //     outputFilePtr() << endl;
        // // }

        // // // write laplacian of snapshots matrix
        // // dataFile = mesh.time().path()/"SVD"/"subDLapFieldValue" + name(zoneI);
        // // outputFilePtr.reset(new OFstream(dataFile));
        // // for (label row = 0; row < subDLapFieldValue.m(); ++row)
        // // {
        // //     for (label column = 0; column < subDLapFieldValue.n(); ++column)
        // //     {
        // //         outputFilePtr().width(16);
        // //         outputFilePtr() << subDLapFieldValue[row][column];
        // //     }
        // //     outputFilePtr() << endl;
        // // }

        // // write subDtemporalCoeff matrix
        // RectangularMatrix<scalar> subDtemporalCoeff(subDFieldValue.n(), modesNum);
        // subDtemporalCoeff = subDFieldValue.T() * subDModesM;
        // dataFile = mesh.time().path()/"SVD"/"temporalCoeff" + name(zoneI);
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

        // // write eigen value matrix
        // dataFile = mesh.time().path()/"SVD"/"eigenVMatrix" + name(zoneI);
        // outputFilePtr.reset(new OFstream(dataFile));
        // for (label row = 0; row < subDfieldValueSVD.S().size(); ++row)
        // {
        //     outputFilePtr() << subDfieldValueSVD.S()[row];
        //     outputFilePtr() << endl;
        // }
