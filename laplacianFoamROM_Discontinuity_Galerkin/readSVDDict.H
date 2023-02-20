
    const word dictName("svdDict");

    IOdictionary svdDict
    (
        IOobject
        (
            dictName, // name of the file
            mesh.time().system(), // path to where the file is
            mesh, // reference to the mesh needed by the constructor
            IOobject::MUST_READ, // indicate that reading this dictionary is compulsory
            IOobject::NO_WRITE
        )

    );
    
    label modesNum (svdDict.lookupOrDefault<label>("modeExtract", mesh.cellZones().size()));
    fileName refCaseName(svdDict.lookupOrDefault<fileName>("modeExtract", "Block1"));
    