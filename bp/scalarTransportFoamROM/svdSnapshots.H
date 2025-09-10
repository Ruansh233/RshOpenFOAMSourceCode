    // create the SVD folder if it is not exist.
    if(!isDir(runTime.globalPath()/"SVD"))
        mkDir(runTime.globalPath()/"SVD");
    
    // ===========================================================
    // ------ some modification to parallel cases are need -------
    // ===========================================================

    // -----------------------------------------------------------
    // --------------------- SVD of T ----------------------------
    // -----------------------------------------------------------
    // SVD for scalar field, T
    // calculate and write matrix element value of diffusion term
    if(mesh.cellZones().size() == 0)
    {
        FatalErrorInFunction
           << "no cellzones in your mesh, please check!\n";
    }

    label elementCellSize (mesh.cellZones()[0].size());
    label snapshotsNum (mesh.cellZones().size());
    RectangularMatrix<scalar> TFieldValue(elementCellSize, snapshotsNum);

    forAll(mesh.cellZones(), zoneI)
    {
        forAll(mesh.cellZones()[zoneI], cellI)
        {
            label cell = mesh.cellZones()[zoneI][cellI];
            TFieldValue(cellI, zoneI) = T[cell];
        }
    }

    SVD TFieldValueSVD(TFieldValue); 
    RectangularMatrix<scalar> TCoefficientsM(TFieldValue.T() * TFieldValueSVD.U());
    RectangularMatrix<scalar> TInvCoefficientsM(SVDinv(TCoefficientsM.T()));


    // -----------------------------------------------------------
    // ---------- ptrlist group for each boundary ----------------
    // -----------------------------------------------------------
    // list of boundary type wordRe
    // this list can read from dictionary
    List<word> boundaryWordRex ({"_in", "_out", "_wall1", 
                                "_wall2", "_wall3", "_wall4"});
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

    Info << "matchPatchID: \n" <<  matchPatchID << endl;
   
    // -----------------------------------------------------------
    // ---------- assign modes value for boundaries --------------
    // -----------------------------------------------------------
    // modes value for each boundary of T, scalar field
    // PtrList can be instead of list, i.e. List<RectangularMatrix<scalar>>
    PtrList<RectangularMatrix<scalar>> TBoundarySnapshotsMList;
    PtrList<RectangularMatrix<scalar>> TBoundaryModesMList;
    RectangularMatrix<scalar> TPatchSnapshots;

    forAll(matchPatchID, patchTypeI)
    {
        TPatchSnapshots.resize(mesh.boundary()[matchPatchID[patchTypeI][0]].size(), snapshotsNum);
        
        forAll(matchPatchID[patchTypeI], patchI)
        {
            forAll(mesh.boundary()[matchPatchID[patchTypeI][patchI]], faceI)
            {
                TPatchSnapshots(faceI, patchI) = T.boundaryField()[matchPatchID[patchTypeI][patchI]][faceI];
            }
        }

        TBoundarySnapshotsMList.append(TPatchSnapshots.clone());
        TBoundaryModesMList.append(RectangularMatrix<scalar> (TPatchSnapshots * TInvCoefficientsM).clone());
    }

    forAll(TBoundarySnapshotsMList, bounI)
    {
        // write snapshots Matrix
        dataFile = runTime.globalPath()/"SVD"/"TBoundarySnapshots"+name(bounI);
        writeMatrix(TBoundarySnapshotsMList[bounI], dataFile);

        // write modes Matrix
        dataFile = runTime.globalPath()/"SVD"/"TBoundaryModes"+name(bounI);
        writeMatrix(TBoundaryModesMList[bounI], dataFile);
    }


    // -----------------------------------------------------------
    // --------- create new object for reference cases -----------
    // -----------------------------------------------------------
    // reference case name
    fileName refCaseName(svdDict.getWord("refCaseName"));
    // time object for reference cases
    Foam::Time runTimeTest
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        // refCaseName is ambigous if it type is word, so the fileName used
        refCaseName,
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
    List<fileName> TModeNames (modesNum);
    forAll(TModeNames, nameI)
    {
        TModeNames[nameI] = "TMode" + name(nameI + 1);
    }

    // -----------------------------------------------------------
    // --------- assign modes value to reference mesh ------------
    // -----------------------------------------------------------
    // ref mesh patch name match ID
    // --- find the patch ID of reference mesh for each patch type
    List<label> refMatchPatchID;

    forAll(boundaryWordRex, RexI)
    {
        matchPatch = ".*" + boundaryWordRex[RexI];
        matchPatch.compile();

        forAll(refElementMesh.boundaryMesh(), patchI)
        {
            if(matchPatch.match(refElementMesh.boundary()[patchI].name()))
            {
                refMatchPatchID.append(patchI);
            }
        }
    }

    Info << "refMatchPatchID: \n" << refMatchPatchID << endl;

    // set mode field value for pressure
    forAll(TModeNames, No_)
    {
        volScalarField TFieldValueMode
        (
            IOobject
            (
                TModeNames[No_],
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
            TFieldValueMode[cellI] = TFieldValueSVD.U()(cellI, No_);
        }

        // set patch value
        forAll(refMatchPatchID, patchI)
        {
            forAll(refElementMesh.boundary()[refMatchPatchID[patchI]], faceI)
            {
                TFieldValueMode.boundaryFieldRef()[refMatchPatchID[patchI]][faceI] =
                    TBoundaryModesMList[patchI](faceI, No_);
            }
        }

        TFieldValueMode.write();
    }


    // write snapshots Matrix
    dataFile = runTime.globalPath()/"SVD"/"TSnapshotsM";
    writeMatrix(TFieldValue, dataFile);

    // write modes Matrix
    dataFile = runTime.globalPath()/"SVD"/"TModesM";
    writeMatrix(TFieldValueSVD.U(), dataFile);

    // write eigen value Matrix
    dataFile = runTime.globalPath()/"SVD"/"TEigenValueM";
    writeMatrix(TFieldValueSVD.S(), dataFile);

    // write eigen value Matrix
    dataFile = runTime.globalPath()/"SVD"/"TCoefficientsM";
    writeMatrix(TCoefficientsM, dataFile);


