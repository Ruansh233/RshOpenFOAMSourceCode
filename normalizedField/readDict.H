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

    const scalar rodWireR(0.5*(rodD+wireD));
    const scalar inclAng(Foam::atan((2*constant::mathematical::pi*rodWireR)/wireH));
    