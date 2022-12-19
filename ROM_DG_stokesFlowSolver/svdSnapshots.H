    // create the SVD folder if it is not exist.
    if(!isDir(mesh.time().path()/"SVD"))
        mkDir(mesh.time().path()/"SVD");
    
    // -----------------------------------------------------------
    // --------------------- SVD of p and U ----------------------
    // -----------------------------------------------------------
    // SVD for scalar field, p
    // calculate and write matrix element value of diffusion term
    label elementCellSize (mesh.cellZones()[0].size());
    label snapshotsNum (mesh.cellZones().size());
    RectangularMatrix<scalar> pFieldValue(elementCellSize, snapshotsNum);

    forAll(mesh.cellZones(), zoneI)
    {
        forAll(mesh.cellZones()[zoneI], cellI)
        {
            label cell = mesh.cellZones()[zoneI][cellI];
            pFieldValue(cellI, zoneI) = p[cell];
        }
    }

    SVD pFieldValueSVD(pFieldValue); 
    RectangularMatrix<scalar> pCoefficientsM(pFieldValue.T() * pFieldValueSVD.U());
    RectangularMatrix<scalar> pInvCoefficientsM(SVDinv(pCoefficientsM.T()));


    // SVD for vector field, U
    // calculate and write matrix element value of diffusion term
    RectangularMatrix<scalar> uFieldValue(3 * elementCellSize, snapshotsNum);

    forAll(mesh.cellZones(), zoneI)
    {
        forAll(mesh.cellZones()[zoneI], cellI)
        {
            label cell = mesh.cellZones()[zoneI][cellI];
            uFieldValue(cellI,                       zoneI) = U[cell].x();
            uFieldValue(cellI + elementCellSize,     zoneI) = U[cell].y();
            uFieldValue(cellI + 2 * elementCellSize, zoneI) = U[cell].z();
        }
    }

    SVD uFieldValueSVD(uFieldValue); 
    RectangularMatrix<scalar> uCoefficientsM(uFieldValue.T() * uFieldValueSVD.U());
    RectangularMatrix<scalar> uInvCoefficientsM(SVDinv(uCoefficientsM.T()));

    // -----------------------------------------------------------
    // ---------- ptrlist group for each boundary ----------------
    // -----------------------------------------------------------
    // list of boundary type wordRe
    // this list can read from dictionary
    List<word> boundaryWordRex ({".*_in", ".*_out", ".*_wall1", 
                                ".*_wall2", ".*_wall3", ".*_wall4"});
    PtrList<List<label>> matchPatchID;
    List<label> tempList;
    wordRe matchPatch;

    
    // the ptr list for each boundary type
    forAll(boundaryWordRex, RexI)
    {
        forAll(mesh.cellZones(), zoneI)
        {
            // matchPatch = boundaryWordRex[RexI];
            // Rsh, 2022-12-19, the outlet of whole domain would be put in front of all interfaces
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
   
    // -----------------------------------------------------------
    // ---------- assign modes value for boundaries --------------
    // -----------------------------------------------------------
    // modes value for each boundary of p, scalar field
    // PtrList can be instead of list, i.e. List<RectangularMatrix<scalar>>
    PtrList<RectangularMatrix<scalar>> pBoundarySnapshotsMList;
    PtrList<RectangularMatrix<scalar>> pBoundaryModesMList;
    RectangularMatrix<scalar> pPatchSnapshots;

    forAll(matchPatchID, patchTypeI)
    {
        pPatchSnapshots.resize(mesh.boundary()[matchPatchID[patchTypeI][0]].size(), snapshotsNum);
        
        forAll(matchPatchID[patchTypeI], patchI)
        {
            forAll(mesh.boundary()[matchPatchID[patchTypeI][patchI]], faceI)
            {
                pPatchSnapshots(faceI, patchI) = p.boundaryField()[matchPatchID[patchTypeI][patchI]][faceI];
            }
        }

        pBoundarySnapshotsMList.append(pPatchSnapshots.clone());
        pBoundaryModesMList.append(RectangularMatrix<scalar> (pPatchSnapshots * pInvCoefficientsM).clone());
    }

    forAll(pBoundarySnapshotsMList, listI)
    {
        dataFile = mesh.time().path()/"SVD"/"pBoundarySnapshots" + name(listI);
        writeMatrix(pBoundarySnapshotsMList[listI], dataFile);

        dataFile = mesh.time().path()/"SVD"/"pBoundaryModes" + name(listI);
        writeMatrix(pBoundaryModesMList[listI], dataFile);
    }


    // modes value for each boundary of U, vector field
    PtrList<RectangularMatrix<scalar>> uBoundarySnapshotsMList;
    PtrList<RectangularMatrix<scalar>> uBoundaryModesMList;
    RectangularMatrix<scalar> uPatchSnapshots;

    forAll(matchPatchID, patchTypeI)
    {
        label patchTypeISize = mesh.boundary()[matchPatchID[patchTypeI][0]].size();
        uPatchSnapshots.resize(3 * patchTypeISize, snapshotsNum);
        
        forAll(matchPatchID[patchTypeI], patchI)
        {
            forAll(mesh.boundary()[matchPatchID[patchTypeI][patchI]], faceI)
            {
                uPatchSnapshots(faceI,                      patchI) = U.boundaryField()[matchPatchID[patchTypeI][patchI]][faceI].x();
                uPatchSnapshots(faceI + patchTypeISize,     patchI) = U.boundaryField()[matchPatchID[patchTypeI][patchI]][faceI].y();
                uPatchSnapshots(faceI + 2 * patchTypeISize, patchI) = U.boundaryField()[matchPatchID[patchTypeI][patchI]][faceI].z();
            }
        }

        uBoundarySnapshotsMList.append(uPatchSnapshots.clone());
        uBoundaryModesMList.append(RectangularMatrix<scalar> (uPatchSnapshots * uInvCoefficientsM).clone());
    }

    forAll(uBoundarySnapshotsMList, listI)
    {
        dataFile = mesh.time().path()/"SVD"/"uBoundarySnapshots" + name(listI);
        writeMatrix(uBoundarySnapshotsMList[listI], dataFile);

        dataFile = mesh.time().path()/"SVD"/"uBoundaryModes" + name(listI);
        writeMatrix(uBoundaryModesMList[listI], dataFile);
    }

    // -----------------------------------------------------------
    // --------- create new object for reference cases -----------
    // -----------------------------------------------------------
    // reference case name
    word refCaseName(svdDict.getWord("refCaseName"));
    // time object for reference cases
    Foam::Time runTimeTest
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        // why refCaseName is ambigous
        "Block1",
        "system",
        "constant"
    );

    // create new mesh object for reference cases
    fvMesh  refElementMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            args.rootPath()/refCaseName/"constant",
            runTimeTest,
            IOobject::MUST_READ
        ),
        false
    );

    // -----------------------------------------------------------
    // ------------------ name for mode file ---------------------
    // -----------------------------------------------------------
    List<fileName> pModeNames (modesNum);
    List<fileName> uModeNames (modesNum);
    forAll(pModeNames, nameI)
    {
        pModeNames[nameI] = "pMode" + name(nameI + 1);
        uModeNames[nameI] = "uMode" + name(nameI + 1);
    }

    // -----------------------------------------------------------
    // --------- assign modes value to reference mesh ------------
    // -----------------------------------------------------------
    // ref mesh patch name match ID
    // --- find the patch ID of reference mesh for each patch type
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

    // set mode field value for pressure
    forAll(pModeNames, No_)
    {
        volScalarField pFieldValueMode
        (
            IOobject
            (
                pModeNames[No_],
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
            pFieldValueMode[cellI] = pFieldValueSVD.U()(cellI, No_);
        }

        // set patch value
        forAll(refMatchPatchID, patchI)
        {
            forAll(refElementMesh.boundary()[refMatchPatchID[patchI]], faceI)
            {
                pFieldValueMode.boundaryFieldRef()[refMatchPatchID[patchI]][faceI] =
                    pBoundaryModesMList[patchI](faceI, No_);
            }
        }

        pFieldValueMode.write();
    }

    // set mode field value for velocity
    forAll(pModeNames, No_)
    {
        volVectorField uFieldValueMode
        (
            IOobject
            (
                uModeNames[No_],
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
            uFieldValueMode[cellI] = Vector<scalar> (uFieldValueSVD.U()(cellI,                     No_),
                                                     uFieldValueSVD.U()(cellI + elementCellSize,   No_),
                                                     uFieldValueSVD.U()(cellI + 2*elementCellSize, No_));
        }

        // set patch value
        forAll(refMatchPatchID, patchI)
        {
            label patchTypeISize = refElementMesh.boundary()[refMatchPatchID[patchI]].size();
            forAll(refElementMesh.boundary()[refMatchPatchID[patchI]], faceI)
            {
                uFieldValueMode.boundaryFieldRef()[refMatchPatchID[patchI]][faceI] = Vector<scalar> (
                    uBoundaryModesMList[patchI](faceI,                      No_),
                    uBoundaryModesMList[patchI](faceI + patchTypeISize,     No_),
                    uBoundaryModesMList[patchI](faceI + 2 * patchTypeISize, No_));
            }
        }

        uFieldValueMode.write();
    }  

    // write snapshots Matrix
    dataFile = mesh.time().path()/"SVD"/"pSnapshotsM";
    writeMatrix(pFieldValue, dataFile);
    dataFile = mesh.time().path()/"SVD"/"uSnapshotsM";
    writeMatrix(uFieldValue, dataFile);

    // write modes Matrix
    dataFile = mesh.time().path()/"SVD"/"pModesM";
    writeMatrix(pFieldValueSVD.U(), dataFile);
    dataFile = mesh.time().path()/"SVD"/"uModesM";
    writeMatrix(pFieldValueSVD.U(), dataFile);

    // write eigen value Matrix
    dataFile = mesh.time().path()/"SVD"/"pEigenValueM";
    writeMatrix(pFieldValueSVD.S(), dataFile);
    dataFile = mesh.time().path()/"SVD"/"uEigenValueM";
    writeMatrix(pFieldValueSVD.S(), dataFile);

    // write eigen value Matrix
    dataFile = mesh.time().path()/"SVD"/"pCoefficientsM";
    writeMatrix(pCoefficientsM, dataFile);
    dataFile = mesh.time().path()/"SVD"/"uCoefficientsM";
    writeMatrix(uCoefficientsM, dataFile);


