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
    DTT1 * T + DTT2
);

forAll(DTT.boundaryField(), patchI)
{
    forAll(DTT.boundaryField()[patchI], faceI)
    {
        DTT.boundaryFieldRef()[patchI][faceI] = DTT1 * T.boundaryFieldRef()[patchI][faceI] + DTT2;
    }
}