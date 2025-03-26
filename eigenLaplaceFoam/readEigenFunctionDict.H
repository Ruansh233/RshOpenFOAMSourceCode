    const word dictName("eigenFunctionDict");

    IOdictionary eigenFunctionDict
    (
        IOobject
        (
            dictName, // name of the file
            mesh.time().system(), // path to where the file is
            mesh, // reference to the mesh needed by the constructor
            IOobject::READ_IF_PRESENT, // indicate that reading this dictionary is compulsory
            IOobject::NO_WRITE
        )

    );

    // Number of eigenfunctions to compute
    label nEigen = eigenFunctionDict.lookupOrDefault<label>("nEigen", 5);
    // Number of iterations
    label maxIter = eigenFunctionDict.lookupOrDefault<label>("maxIter", 1000);
    // Tolerance
    scalar tol = eigenFunctionDict.lookupOrDefault<scalar>("tol", 1e-6);
    