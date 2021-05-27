void MMC_tPhi(double Cxy00[2], double dx, double dy, double cellpos_x[], double cellpos_z[], int N, double tphi[], double p, double hi)
{
    double *x = new double[N];
    double *y = new double[N];
    int i;
    double st = std::sin(0);
    double ct = std::cos(0);
    for (i = 0; i < N; i++)
    {
        x[i] = ct * (cellpos_x[i] - Cxy00[0]) + st * (cellpos_z[i] - Cxy00[1]);
        y[i] = -st * (cellpos_x[i] - Cxy00[0]) + ct * (cellpos_z[i] - Cxy00[1]);
        tphi[i] = MMC_Heaviside(1 - std::pow(x[i] / dx, p) - std::pow(y[i] / dy, p), hi);
    }
    delete x;
    delete y;
}
