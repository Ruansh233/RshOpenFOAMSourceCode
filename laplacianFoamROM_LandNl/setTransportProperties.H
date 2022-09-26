forAll(mesh.C(), cellI)
{
    DT[cellI] = DTT1 * T[cellI] + DTT2;
}
// Info << "DTT1: " << DTT1 * max(T).value() << ", " << "DTT2: " << DTT2 << ".\n";
// Info << "Nonlinear ratio: " << DTT1 * max(T).value() / DTT2 << ". \n";
