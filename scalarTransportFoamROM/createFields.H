Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Info<< "Reading field U\n" << endl;

// volVectorField U
// (
//     IOobject
//     (
//         "U",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );

Info<< "Creating field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimVelocity
);

wallDist wallDistance_(mesh);

forAll(U, cellI)
{
    U[cellI] = vector(0,0,wallDistance_.y()[cellI]);
}

forAll(U.boundaryField(), patchI)
{
    U.boundaryFieldRef().set(patchI, 
        fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], U));
}

U.write();

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


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT("DT", dimViscosity, transportProperties);

#include "createPhi.H"

#include "createFvOptions.H"
