    word subdomainName (customDict.getWord("subdomainName"));
    word subdomainZoneList (customDict.getWord("subdomainZoneList"));
    word dataPath (customDict.getWord("dataPath"));
    
    List<word> scalarZoneFieldName (customDict.lookup("scalarZoneFieldName"));
    List<word> vectorZoneFieldName (customDict.lookup("vectorZoneFieldName"));

    // List<List<label>> subdomains;
    
    // output the subdomains List
    IOList<List<label>> subdomainsIO
    (
        IOobject
        (
            subdomainZoneList,
            runTime.caseConstant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
