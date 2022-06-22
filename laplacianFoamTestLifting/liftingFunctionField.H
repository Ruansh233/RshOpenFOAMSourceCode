// Field<scalar> liftingFunctionField (mesh.C().size());
wordList liftPatchTypes_ ({"zeroGradient", "zeroGradient", "zeroGradient", "zeroGradient", 
                            "fixedValue", "fixedValue"});

// create field 
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
    dimless,
    liftPatchTypes_
);

// assign cell value
forAll(mesh.C(), cellI)
{
    liftingFunction[cellI] = 273 + (mesh.C()[cellI].component(vector::Z) + 0.2) / 0.4 * (573 - 273); 
}

// assign boundary value
forAll(liftingFunction.boundaryField(), patchI)
{
    forAll(liftingFunction.boundaryField()[patchI], faceI)
    {
        liftingFunction.boundaryFieldRef()[patchI][faceI] = 273 + 
                        (mesh.boundary()[patchI].Cf()[faceI].component(vector::Z) + 0.2) / 0.4 * 
                        (573 - 273);
    }
}

liftingFunction.write();

volVectorField liftGrad (fvc::grad(liftingFunction));
liftGrad.write();
