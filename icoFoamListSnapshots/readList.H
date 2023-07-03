    // read cell List
    const label subdomainNum = 20;
    List<List<label>> globalcellIndexList(subdomainNum);

    autoPtr<IFstream> inputFilePtr;
    inputFilePtr.reset(new IFstream(runTime.caseConstant()/"globalcellIndexList"));
    inputFilePtr() >> globalcellIndexList;

    // read face List
    const label patchNum = 6;
    PtrList<List<List<label>>> globalfaceIndexListPtr;
    List<List<label>> globalfaceIndexList;
    for(label patchI=0; patchI<patchNum; ++patchI)
    {
        inputFilePtr.reset(new IFstream(runTime.caseConstant()/"globalfaceIndexList_"+name(patchI)));
        inputFilePtr() >> globalfaceIndexList;
        
        globalfaceIndexListPtr.append( globalfaceIndexList.clone());
    }
            
