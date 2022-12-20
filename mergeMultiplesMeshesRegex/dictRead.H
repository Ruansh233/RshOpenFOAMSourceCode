    // Get access to a custom dictionary
    // dictionary customDict;
    // const word dictName("mergeMultipleMeshesDict");

    // // Create and input-output object - this holds the path to the dict and its name
    // IOdictionary customDict
    // (
    //     IOobject
    //     (
    //         dictName, // name of the file
    //         mesh.time().system(), // path to where the file is
    //         mesh, // reference to the mesh needed by the constructor
    //         IOobject::MUST_READ // indicate that reading this dictionary is compulsory
    //     )
    // );

    // Rsh, 22.11.22
    // --- add dict option
    const word dictName("mergeMultipleMeshesDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;
    IOdictionary customDict(dictIO);

    fileName masterCase(customDict.lookup("masterCase"));
    List<fileName> addCaseListTemp(customDict.lookup("addCase"));    
    List<fileName> addCaseList;
    List<fileName> foldersList(readDir(mesh.time().rootPath(), fileName::DIRECTORY));
    wordRe fileNameRex;

    // Rsh, 2022-12-17
    // --- add sort of the folder name list
    // --- since sortableList in openfoam doesn't work, sort of std::list is employed
    std::list<std::string> testsortList;
    forAll(foldersList, listI)
    {
        testsortList.push_back(foldersList[listI]);
    }
    testsortList.sort();
    label indexI = 0;
    for (std::string n : testsortList)
    {
        foldersList[indexI] = n;
        ++indexI;
    }        

    forAll(addCaseListTemp, listI)
    {
        fileNameRex = addCaseListTemp[listI];
        fileNameRex.compile();
        
        if(fileNameRex.isPattern())
        { 
            forAll(foldersList, folder)
            {
                if(fileNameRex.match(foldersList[folder]))
                    addCaseList.append(mesh.time().rootPath() + "/" + foldersList[folder]);
            }
        }
        else
        {
            addCaseList.append(addCaseListTemp[listI]);
        }
    }

