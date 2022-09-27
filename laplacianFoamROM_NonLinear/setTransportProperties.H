scalar DTT (transportProperties.getScalar("DTT"));

volScalarField DTTT 
(
    IOobject
    (
        "DTTT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimViscosity,
    DTT * T
);

forAll(DTTT.boundaryField(), patchI)
{
    forAll(DTTT.boundaryField()[patchI], faceI)
    {
        DTTT.boundaryFieldRef()[patchI][faceI] = DTT * T.boundaryFieldRef()[patchI][faceI];
    }
}

// forAll(mesh.C(), cellI)
// {
//     DT[cellI] = DTT * T[cellI];
// }
