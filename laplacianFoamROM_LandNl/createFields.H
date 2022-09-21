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


Info<< "Reading diffusivity DT\n" << endl;

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity, Zero)
);

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

if (!DT.headerOk())
{
    DT = dimensionedScalar("DT", dimViscosity, transportProperties);
}

scalar DTT1 (transportProperties.getScalar("DTT1"));
scalar DTT2 (transportProperties.getScalar("DTT2"));

Info << "DTT1: " << DTT1 << nl
     << "DTT2: " << DTT2 << nl;

#include "createFvOptions.H"
