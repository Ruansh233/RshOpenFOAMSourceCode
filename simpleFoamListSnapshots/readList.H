    // read cell List
    List<List<label>> cellIndexList;
    IOList<List<label>> cellIndexListIO
    (
        IOobject
        (
            "globalcellIndexList",
            runTime.caseConstant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        cellIndexList
    );
    

    // read face List
    const label patchNum = 6;
    PtrList<List<List<label>>> faceIndexListPtr(patchNum);
    List<List<label>> faceIndexList;
    for(label patchI=0; patchI<patchNum; ++patchI)
    {
        IOList<List<label>> faceIndexListIO
        (
            IOobject
            (
                "globalfaceIndexList_"+name(patchI),
                runTime.caseConstant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            faceIndexList
        );
        
        faceIndexListPtr[patchI] = faceIndexList.clone();
    }
    
    Info<< "cellIndexList " << cellIndexList << endl;

    Info<< "faceIndexListPtr " << faceIndexListPtr[0] << endl;
            
