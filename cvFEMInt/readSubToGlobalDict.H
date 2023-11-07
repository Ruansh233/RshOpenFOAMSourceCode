    const word dictName("subToGlobalDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;
    IOdictionary subToGlobalDict(dictIO);

    fileName refCaseName(subToGlobalDict.getWord("refCase"));
    Info<< "Reference case: " << refCaseName << endl;
    
    List<word> scalarFields(subToGlobalDict.lookup("scalarFields"));
    List<word> vectorFields(subToGlobalDict.lookup("vectorFields"));
    List<word> tensorFields(subToGlobalDict.lookup("tensorFields"));
    