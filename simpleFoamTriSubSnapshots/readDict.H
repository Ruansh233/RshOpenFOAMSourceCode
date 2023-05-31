    List<List<label>> subdomains(customDict.getLabel("subdomainNum"));
    List<List<word>> subdomainsZoneName(customDict.getLabel("subdomainNum"));
    List<word> subdomainDict;

    forAll(subdomains, domainI)
    {
        subdomainDict = customDict.subDict("subdomains").lookup("tri_"+name(domainI+1));

        forAll(subdomainDict, subDI)
        {
            label zoneID = mesh.cellZones().findZoneID(subdomainDict[subDI]);
            subdomains[domainI].append(zoneID);
            subdomainsZoneName[domainI].append(subdomainDict[subDI]);
        }
    }
    
    // output the subdomains List
    IOList<List<label>> subdomainsIO
    (
        IOobject
        (
            "subdomains",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subdomains
    );
    subdomainsIO.write();  


    IOList<List<word>> subdomainsZoneNameIO
    (
        IOobject
        (
            "subdomainsZoneName",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subdomainsZoneName
    );
    subdomainsZoneNameIO.write();  