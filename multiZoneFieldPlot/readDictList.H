    word subdomainName (customDict.getWord("subdomainName"));
    List<word> subdomainZonesList (customDict.lookup("subdomainZonesList"));
    List<word> dataPathList (customDict.lookup("dataPathList"));
    
    List<word> scalarZoneFieldName (customDict.lookup("scalarZoneFieldName"));
    List<word> vectorZoneFieldName (customDict.lookup("vectorZoneFieldName"));

    // List<List<label>> subdomains;

    Info<< "The input data is shown following, " << endl
        << "subdomainZonesList: " << subdomainZonesList << endl
        << "dataPathList: " << dataPathList << endl
        << "scalarZoneFieldName: " << scalarZoneFieldName << endl
        << "vectorZoneFieldName: " << vectorZoneFieldName << endl;
