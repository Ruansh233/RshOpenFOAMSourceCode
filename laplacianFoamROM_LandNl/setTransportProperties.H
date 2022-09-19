scalar DTT (transportProperties.getScalar("DTT"));
scalar referenceT (transportProperties.getScalar("referenceT"));

forAll(mesh.C(), cellI)
{
    DT[cellI] = DTT * (T[cellI] - referenceT);
    // DT[cellI] = DTT* T[cellI];
}
