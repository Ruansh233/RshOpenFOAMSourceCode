    const word dictName("MomentumSourceDict");

    // Create and input-output object - this holds the path to the dict and its name
    IOdictionary customDict
    (
        IOobject
        (
            dictName, // name of the file
            mesh.time().system(), // path to where the file is
            mesh, // reference to the mesh needed by the constructor
            IOobject::MUST_READ // indicate that reading this dictionary is compulsory
        )
    );

    const scalar rodD(customDict.getScalar("rodD"));
    const scalar wireD(customDict.getScalar("wireD"));
    const scalar wireH(customDict.getScalar("wireH"));
    const List<vector> rodCentroids (customDict.lookup("rodCentroids"));

    const scalar refHight(customDict.getOrDefault<scalar>("refHight", 0.0));
    const scalar refAngle(customDict.getOrDefault<scalar>("refAngle", 0.0));
    const word cellZoneName;
    customDict.readIfPresent("cellZoneName", cellZoneName);
    