    const word dictName("DGdict");

    const IOdictionary DGdict
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
    const label modesNum (DGdict.getLabel("modesNum"));
    const label elementNum (DGdict.getLabel("elementNum"));
    const label iterations (DGdict.getLabel("iterations"));

    // Discontinuity galerkin method parameters, add to dict reading later
    scalar epsilonPara (DGdict.getScalar("epsilonPara"));
    scalar xigema0 (DGdict.getScalar("xigema0"));
    scalar xigema1 (DGdict.getScalar("xigema1"));
    scalar xigema2 (DGdict.getScalar("xigema2"));
    scalar xigema3 (DGdict.getScalar("xigema3"));
    scalar xigema4 (DGdict.getScalar("xigema4"));
    scalar xigema5 (DGdict.getScalar("xigema5"));
    // scalar beta0 (DGdict.getScalar("beta0"));
    // scalar beta1 (DGdict.getScalar("beta1"));

    // iteration weighting factor
    const scalar weightW (DGdict.getScalar("weightW"));

    // phiscal parameter    
    const scalar heatConductivity (DGdict.getScalar("heatConductivity"));
    // scalar Tin (DGdict.getScalar("Tin"));
    // scalar Tout (DGdict.getScalar("Tout"));
    // scalar qn (DGdict.getScalar("qn"));
    const vector Uin (DGdict.lookup("Uin"));

    // The following is same to the lookup function
    // vector UinTest (DGdict.get<vector>("Uin"));

    


