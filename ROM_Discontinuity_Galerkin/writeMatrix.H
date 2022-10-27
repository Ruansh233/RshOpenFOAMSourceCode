template <class Form, class Type>
void writeMatrix(Matrix< Form, Type > matrixProc,
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

template <class Form, class Type>
void writeMatrix(Matrix< Form, Type > matrixProc,
                fileName dataFile)
{
    return writeMatrix(matrixProc, dataFile,
                        matrixProc.m(),
                        matrixProc.n()
    );
}

template <class Type>
void writeMatrix(List<Type> matrixProc,
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

template <class Type>
void writeMatrix(List<Type> matrixProc,
                fileName dataFile)
{
    return writeMatrix(matrixProc, dataFile,
                        matrixProc.size()
    );
}
