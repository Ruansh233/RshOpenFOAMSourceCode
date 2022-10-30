    const word dictName("DGdict");

    IOdictionary DGdict
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

    // system parameter
    label modesNum (DGdict.getLabel("modesNum"));
    label elementNum (DGdict.getLabel("elementNum"));

    // Discontinuity galerkin method parameters, add to dict reading later
    scalar epsilonPara (DGdict.getScalar("epsilonPara"));
    scalar xigema0 (DGdict.getScalar("xigema0"));
    scalar xigema1 (DGdict.getScalar("xigema1"));
    scalar beta0 (DGdict.getScalar("beta0"));
    scalar beta1 (DGdict.getScalar("beta1"));

    // phiscal parameter    
    scalar heatConductivity (DGdict.getScalar("heatConductivity"));