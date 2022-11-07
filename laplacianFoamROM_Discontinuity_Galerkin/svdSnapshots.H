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
    // this list can read from dictionary
    List<word> boundaryWordRex ({".*_in", ".*_out", ".*_wall1", ".*_wall2", ".*_wall3", ".*_wall4"});
    PtrList<List<label>> matchPatchID;
    List<label> tempList;
    wordRe matchPatch;

    
    // the ptr list for each boundary type
    forAll(boundaryWordRex, RexI)
    {
        forAll(mesh.cellZones(), zoneI)
        {
            // matchPatch = boundaryWordRex[RexI];
            // the outlet of whole domain would be put in front of all interfaces
            // --- so name(zoneI + 1) is used
            matchPatch = ".*block" + name(zoneI + 1) + boundaryWordRex[RexI];
            matchPatch.compile();

            forAll(mesh.boundaryMesh(), patchI)
            {
                if(matchPatch.match(mesh.boundary()[patchI].name()))
                {
                    tempList.append(patchI);
                }
            }
        }

        matchPatchID.append(tempList.clone());
        tempList.clear();
    }
   

    // modes value for each boundary
    // PtrList can be instead of list, i.e. List<RectangularMatrix<scalar>>
    PtrList<RectangularMatrix<scalar>> boundarySnapshotsMList;
    PtrList<RectangularMatrix<scalar>> boundaryModesMList;
    RectangularMatrix<scalar> patchSnapshots;

    forAll(matchPatchID, patchTypeI)
    {
        patchSnapshots.resize(mesh.boundary()[matchPatchID[patchTypeI][0]].size(), snapshotsNum);
        
        forAll(matchPatchID[patchTypeI], patchI)
        {
            forAll(mesh.boundary()[matchPatchID[patchTypeI][patchI]], faceI)
            {
                patchSnapshots(faceI, patchI) = T.boundaryField()[matchPatchID[patchTypeI][patchI]][faceI];
            }
        }

        boundarySnapshotsMList.append(patchSnapshots.clone());
        boundaryModesMList.append(RectangularMatrix<scalar> (patchSnapshots * invCoefficientsM).clone());
    }

    forAll(boundarySnapshotsMList, listI)
    {
        dataFile = mesh.time().path()/"SVD"/"boundarySnapshots" + name(listI);
        writeMatrix(boundarySnapshotsMList[listI], dataFile);

        dataFile = mesh.time().path()/"SVD"/"boundaryModes" + name(listI);
        writeMatrix(boundaryModesMList[listI], dataFile);
    }


    // create new time object for other cases
    Foam::Time runTimeTest
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        // args.caseName(),
        "Block1",
        "system",
        "constant"
    );

    // create new mesh object for other cases
    fvMesh  refElementMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            args.rootPath()/"Block1"/"constant",
            runTimeTest,
            IOobject::MUST_READ
        ),
        false
    );


    forAll(modeNames, nameI)
    {
        modeNames[nameI] = "Tmode" + name(nameI + 1);
    }


    // ref mesh patch name match ID
    List<label> refMatchPatchID;

    forAll(boundaryWordRex, RexI)
    {
        matchPatch = boundaryWordRex[RexI];
        matchPatch.compile();

        forAll(refElementMesh.boundaryMesh(), patchI)
        {
            if(matchPatch.match(refElementMesh.boundary()[patchI].name()))
            {
                refMatchPatchID.append(patchI);
            }
        }
    }


    // set mode field value
    forAll(modeNames, No_)
    {
        // create mode field by copying T
        volScalarField fieldValueMode
        (
            IOobject
            (
                modeNames[No_],
                runTimeTest.timeName(),
                refElementMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            refElementMesh,
            dimless
        );

        // set cell value
        forAll(refElementMesh.C(), cellI)
        {
            fieldValueMode[cellI] = fieldValueSVD.U()(cellI, No_);
        }

        // set patch value
        forAll(refMatchPatchID, patchI)
        {
            forAll(refElementMesh.boundary()[refMatchPatchID[patchI]], faceI)
            {
                fieldValueMode.boundaryFieldRef()[refMatchPatchID[patchI]][faceI] =
                    boundaryModesMList[patchI](faceI, No_);
            }
        }

        fieldValueMode.write();

        // create mode field by copying T
        volVectorField fieldValueModegrad
        (
            IOobject
            (
                "grad" + modeNames[No_],
                runTimeTest.timeName(),
                refElementMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(fieldValueMode)
        );
        fieldValueModegrad.write();
    }


    // // grad and laplacian of element fieldValueMode
    // forAll(modeNames, No_)
    // {
    //     // create mode field by copying T
    //     volScalarField fieldValueModeTest
    //     (
    //         IOobject
    //         (
    //             modeNames[No_],
    //             runTimeTest.timeName(),
    //             refElementMesh,
    //             IOobject::READ_IF_PRESENT,
    //             IOobject::AUTO_WRITE
    //         ),
    //         refElementMesh
    //     );

    //     // create mode field by copying T
    //     volVectorField fieldValueModegrad
    //     (
    //         IOobject
    //         (
    //             "grad" + modeNames[No_],
    //             runTimeTest.timeName(),
    //             refElementMesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         fvc::grad(fieldValueModeTest)
    //     );
    //     fieldValueModegrad.write();

    //     // create mode field by copying T
    //     volScalarField fieldValueModelap
    //     (
    //         IOobject
    //         (
    //             "lap" + modeNames[No_],
    //             runTimeTest.timeName(),
    //             refElementMesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         fvc::laplacian(fieldValueModeTest)
    //     );
    //     fieldValueModelap.write();
    // }
    

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


