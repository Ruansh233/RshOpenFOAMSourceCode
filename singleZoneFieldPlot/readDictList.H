    word subdomainName (customDict.getWord("subdomainName"));
    word cellZonesType (customDict.getWord("cellZonesType"));
    word dataPath (customDict.getWord("dataPath"));
    List<word> typesName(customDict.lookup("types"));
    
    List<List<word>> scalarZoneFieldName (customDict.lookup("scalarZoneFieldName"));
    List<List<word>> vectorZoneFieldName (customDict.lookup("vectorZoneFieldName"));

    // List<List<label>> subdomains;
    
    // output the subdomains List
    IOList<List<label>> cellZonesTypeIO
    (
        IOobject
        (
            cellZonesType,
            runTime.caseConstant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Info<< scalarZoneFieldName << endl;

    Info<< "The input data is shown following, " << endl
        << "subdomainName: " << subdomainName << endl
        << "dataPath: " << dataPath << endl
        << "scalarZoneFieldName: " << scalarZoneFieldName << endl
        << "vectorZoneFieldName: " << vectorZoneFieldName << endl;
