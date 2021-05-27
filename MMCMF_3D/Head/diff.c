double diff(double x[], double del, double eta, int allcells)
{
    int i;
    double z = 0;
    double *xh = new double[allcells];
    for (i = 0; i < allcells; i++)
    {
        if (x[i] <= eta)
        {
            xh[i] = eta * (Foam::exp(-del * (1 - x[i] / eta)) - (1 - x[i] / eta) * Foam::exp(-del));
        }
        else
        {
            xh[i] = eta + (1 - eta) * (1 - Foam::exp(-del * (x[i] - eta) / (1 - eta)) + (x[i] - eta) * Foam::exp(-del) / (1 - eta));
        }
    }
    for (i = 0; i < allcells; i++)
    {
        z = z + x[i] - xh[i];
    }
    delete xh;
    return {z};
}
