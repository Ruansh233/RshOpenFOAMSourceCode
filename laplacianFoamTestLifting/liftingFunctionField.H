Field<scalar> liftingFunctionField (mesh.C().size());

forAll(mesh.C(), cellI)
{
    liftingFunctionField[cellI] = 273 + (mesh.C()[cellI].component(vector::Z) + 0.2) / 0.4 * (573 - 273); 
}

volScalarField liftingFunction
(
    IOobject
    (
        "liftingFunction",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimTemperature,
    liftingFunctionField
);

forAll(liftingFunction.boundaryField(), patchI)
{
    forAll(liftingFunction.boundaryField()[patchI], faceI)
    {
        liftingFunction.boundaryFieldRef()[patchI][faceI] = 273 + (mesh.boundary()[patchI].Cf()[faceI].component(vector::Z) + 0.2) / 0.4 * (573 - 273);
    }
}

// liftingFunction.write();
