volScalarField DTlapT
(
    IOobject
    (
        "DTlapT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::laplacian(DTTT, T)
);

volScalarField DT_lapT
(
    IOobject
    (
        "DT_lapT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    DTTT * fvc::laplacian(T)
);

volScalarField Div_gradT
(
    IOobject
    (
        "Div_gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(gradT)
);

volScalarField Div_DT_gradT
(
    IOobject
    (
        "Div_DT_gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::div(DTTT * gradT)
);

volScalarField gradDT_gradT
(
    IOobject
    (
        "gradDT_gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(DTTT) & gradT
);

volScalarField gradDT_gradT_DT_lapT
(
    IOobject
    (
        "gradDT_gradT_DT_lapT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    gradDT_gradT + DT_lapT
);