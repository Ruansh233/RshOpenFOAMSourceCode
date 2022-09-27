volScalarField DTT
(
    IOobject
    (
        "DTT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimViscosity,
    DK * T
);

forAll(DTT.boundaryField(), patchI)
{
    forAll(DTT.boundaryField()[patchI], faceI)
    {
        DTT.boundaryFieldRef()[patchI][faceI] = DK * T.boundaryFieldRef()[patchI][faceI];
    }
}
