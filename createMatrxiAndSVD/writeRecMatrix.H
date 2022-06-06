void wirteMatrix(Matrix< Form, Type > &outMatrix, fileName &outputDir)
{
    OFstream outputFile_(outputDir);

    labelPair matrixSizes_ (outMatrix.sizes());
    if (matrixSizes_.second() == 1)
    {
        for (label row = 0; row < outMatrix.size(); ++row)
        {
            outputFile_ << outMatrix[row];
            outputFile_ << endl;
        }
    }
    else
    {
        for (label row = 0; row < outMatrix.m(); ++row)
        {
            for (label column = 0; column < outMatrix.n(); ++column)
            {
                outputFile_.width(16);
                outputFile_ << outMatrix[row][column];
            }
            outputFile_ << endl;
        }
    }
}