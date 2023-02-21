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
    mesh
);

wallDist wallDistance_(mesh);

forAll(U, cellI)
{
    U[cellI] = vector(0,0,wallDistance_.y()[cellI]);
}

