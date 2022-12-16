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

    List<scalar> testList1 ({1.0, 2.0, 3.0});

    SortableList<scalar> testList2(List<scalar> (testList1));
    // testList2.clean ();  

    Info << "\nfolders: " << testList2 << endl;
