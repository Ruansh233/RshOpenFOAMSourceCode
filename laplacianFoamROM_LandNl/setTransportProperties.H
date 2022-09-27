dimensionedScalar DT1 = dimensionedScalar("DT1", dimensionSet(0, 2, -1, -1, 0, 0, 0), DTT1);
dimensionedScalar DT2 = dimensionedScalar("DT2", dimViscosity, DTT2);

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
    DT1 * T + DT2
);

forAll(DTT.boundaryField(), patchI)
{
    forAll(DTT.boundaryField()[patchI], faceI)
    {
        DTT.boundaryFieldRef()[patchI][faceI] = DT1.value() * T.boundaryFieldRef()[patchI][faceI] + DT2.value();
    }
}