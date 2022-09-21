forAll(mesh.C(), cellI)
{
    DT[cellI] = DTT1 * T[cellI] + DTT2;
}
