volScalarField yPlusTemp
(
    IOobject
    (
        "yPlusTemp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("yPlusTemp", dimless, 0.0)
);

volScalarField yPlus
(
    IOobject
    (
        "yPlus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("yPlus", dimless, 0.0)
);

volVectorField uPlus
(
    IOobject
    (
        "uPlus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("uPlus", dimless, vector::zero)
);

volScalarField uTau
(
    IOobject
    (
        "uTau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("uTau", dimVelocity, SMALL)
);

Info<< "Reading field U\n" << endl;
volVectorField U
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

#include "createPhi.H"


// label pRefCell = 0;
// scalar pRefValue = 0.0;
// setRefCell(p, simple.dict(), pRefCell, pRefValue);
// mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

volScalarField::Boundary d = nearWallDist(mesh).y();
