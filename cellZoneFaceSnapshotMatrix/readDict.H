    List<word> typesName(customDict.lookup("types"));
    label casesNumber(customDict.getLabel("CasesNumber"));
    List<label> CasesID(customDict.lookup("CasesID"));

    word typeSubDictName;
    word caseSubDictName;
    List<wordRe> includeZones;
    List<wordRe> excludeZones;

    List<List<List<label>>> cellZonesType(typesName.size(), List<List<label>> (casesNumber));
    List<List<List<word>>> cellZonesTypeName(typesName.size(), List<List<word>> (casesNumber));

    forAll(typesName, typeI)
    {
        typeSubDictName = typesName[typeI];
        for(label i=0; i<casesNumber; ++i)
        {
            caseSubDictName = "case" + name(CasesID[i]);
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
                                cellZonesType[typeI][i].append(zoneI);
                                cellZonesTypeName[typeI][i].append(mesh.cellZones()[zoneI].name());
                            }                                
                        }
                    }
                }
            }
        }
    }

    
    forAll(cellZonesTypeName, typeI)
    {
        IOList<List<label>> cellZonesTypeIO
        (
            IOobject
            (
                "cellZonesType"+name(typeI),
                runTime.caseConstant(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            cellZonesType[typeI]
        );
        cellZonesTypeIO.write();  


        IOList<List<word>> cellZonesTypeNameIO
        (
            IOobject
            (
                "cellZonesTypeName"+name(typeI),
                runTime.caseConstant(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            cellZonesTypeName[typeI]
        );
        cellZonesTypeNameIO.write();  
    }
    

            
