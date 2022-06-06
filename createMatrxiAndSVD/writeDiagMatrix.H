void wirteMatrix(Matrix< Form, Type > &outMatrix, fileName &outputDir)
{
    OFstream outputFile_(outputDir);
    for (label row = 0; row < outMatrix.size(); ++row)
    {
        outputFile_ << outMatrix[row];
        outputFile_ << endl;
    }
}