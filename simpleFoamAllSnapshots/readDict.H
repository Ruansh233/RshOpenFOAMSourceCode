    List<word> typesName(customDict.lookup("types"));

    word typeSubDictName;
    word caseSubDictName;
    List<wordRe> includeZones;
    List<wordRe> excludeZones;

    if(! customDict.readIfPresent("caseName", caseSubDictName))
    {
        Info<< "The Entry 'caseName' not found in dictionary.\n \
                caseName of current case name is used." << endl;
        caseSubDictName = runTime.caseName();
        Info<< "caseName: " << caseSubDictName << endl;
    }

    List<List<label>> cellZonesType(typesName.size());
    List<List<word>> cellZonesTypeName(typesName.size());

    forAll(typesName, typeI)
    {
        typeSubDictName = typesName[typeI];

        includeZones = customDict.subDict(typeSubDictName).subDict(caseSubDictName).lookup("includeZones");
        excludeZones = customDict.subDict(typeSubDictName).subDict(caseSubDictName).lookup("excludeZones");

        forAll(mesh.cellZones(), zoneI)
        {
            forAll(includeZones, includeI)
            {
                wordRe tmpRegex (includeZones[includeI]);
                tmpRegex.compile();

                if(tmpRegex.match(mesh.cellZones()[zoneI].name()))
                // match of include could be done only 1 time
                {
                    label excludeFlag(0);
                    forAll(excludeZones, excludeI)
                    {
                        tmpRegex = excludeZones[excludeI];
                        tmpRegex.compile();

                        if(tmpRegex.match(mesh.cellZones()[zoneI].name()))
                            break;
                        else
                            ++excludeFlag;

                        if(excludeFlag == excludeZones.size())
                        {
                            cellZonesType[typeI].append(zoneI);
                            cellZonesTypeName[typeI].append(mesh.cellZones()[zoneI].name());
                        }                                
                    }
                }
            }
        }
    }

    
    // output the cellZone List
    IOList<List<label>> cellZonesTypeIO
    (
        IOobject
        (
            "cellZonesType",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellZonesType
    );
    cellZonesTypeIO.write();  


    IOList<List<word>> cellZonesTypeNameIO
    (
        IOobject
        (
            "cellZonesTypeName",
            runTime.caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        cellZonesTypeName
    );
    cellZonesTypeNameIO.write();  
    

            
