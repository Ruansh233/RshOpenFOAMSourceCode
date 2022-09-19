scalar DTT (transportProperties.getScalar("DTT"));

forAll(mesh.C(), cellI)
{
    // DT[cellI] = 1.0e-5 * (600 - T[cellI]);
    DT[cellI] = DTT* T[cellI];
}
