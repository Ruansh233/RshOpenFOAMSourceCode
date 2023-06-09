    List<List<label>> subdomains(customDict.getLabel("subdomainNum"));
    List<List<word>> subdomainsZoneName(customDict.getLabel("subdomainNum"));
    List<word> subdomainDict;
    label domainNum (0);

    forAll(subdomains, domainI)
    {
        if(customDict.subDict("subdomains").readIfPresent("sub_"+name(domainI+1), subdomainDict))
        {
            forAll(subdomainDict, subDI)
            {
                label zoneID = mesh.cellZones().findZoneID(subdomainDict[subDI]);
                subdomains[domainNum].append(zoneID);
                subdomainsZoneName[domainNum].append(subdomainDict[subDI]);
            }
            ++ domainNum; 
        }   
    }

    subdomains.resize(domainNum);
    subdomainsZoneName.resize(domainNum);
    
    word subdomainName(customDict.getWord("subdomainName"));
    // output the subdomains List
    IOList<List<label>> subdomainsIO
    (
        IOobject
        (
            "subdomains_"+subdomainName,
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
            "subdomainsZoneName_"+subdomainName,
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        subdomainsZoneName
    );
    subdomainsZoneNameIO.write();  
