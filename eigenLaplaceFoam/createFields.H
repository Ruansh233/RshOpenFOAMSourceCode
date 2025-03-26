Info<< "Reading field U\n" << endl;

volScalarField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField v
(
    IOobject
    (
        "v", 
        runTime.timeName(), 
        mesh, 
        IOobject::NO_READ, 
        IOobject::NO_WRITE
    ),
    U
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    transportProperties
);

#include "createFvOptions.H"
