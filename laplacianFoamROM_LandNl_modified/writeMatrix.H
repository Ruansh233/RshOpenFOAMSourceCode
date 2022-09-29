void writeMatrix(RectangularMatrix<scalar> matrixProc,
                fileName dataFile,
                label rows,
                label columns)
{
    OFstream writeOF(dataFile);
    for (label row = 0; row < rows; ++row)
    {
        for (label column = 0; column < columns; ++column)
        {
            writeOF.width(16);
            writeOF << matrixProc[row][column];
        }
        writeOF << endl;
    }
}

void writeMatrix(RectangularMatrix<scalar> matrixProc,
                fileName dataFile)
{
    return writeMatrix(matrixProc, dataFile,
                        matrixProc.m(),
                        matrixProc.n()
    );
}
