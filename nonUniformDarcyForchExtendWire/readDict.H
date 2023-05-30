    const word dictName("DFmodelDict");

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
    const scalar wireInfluD(customDict.getScalar("wireInfluD"));
    const scalar wireH(customDict.getScalar("wireH"));
    const List<vector> rodCentroids (customDict.lookup("rodCentroids"));

    const scalar rodWireP(0.5*(rodD+wireInfluD));

    const scalar D0(customDict.getScalar("D0"));
    const scalar F0(customDict.getScalar("F0"));

    const scalar xigmax(customDict.getScalar("xigmax"));
    const scalar xigmay(customDict.getScalar("xigmay"));
    