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
            writeOF << matrixProc(row, column);
            writeOF << " ";
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

void writeMatrix(scalarDiagonalMatrix matrixProc,
                fileName dataFile,
                label rows)
{
    OFstream writeOF(dataFile);
    for (label row = 0; row < rows; ++row)
    {
        writeOF << matrixProc[row];
        writeOF << endl;
    }
}

void writeMatrix(scalarDiagonalMatrix matrixProc,
                fileName dataFile)
{
    return writeMatrix(matrixProc, dataFile,
                        matrixProc.size()
    );
}

template <class T>
void writeMatrix(T matrixProc,
                fileName dataFile,
                label rows)
{
    OFstream writeOF(dataFile);
    for (label row = 0; row < rows; ++row)
    {
        writeOF << matrixProc[row];
        writeOF << endl;
    }
}

template <class T>
void writeMatrix(T matrixProc,
                fileName dataFile)
{
    return writeMatrix(matrixProc, dataFile,
                        matrixProc.size()
    );
}
