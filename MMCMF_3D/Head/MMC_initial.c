void MMC_initial(int CW, int CH, double DW, double DH, double Dxx[], double xy00[], double femxy00[])
{
    int i;
    double x_int = 1.0 / CW;
    double y_int = 1.0 / CH;
    double *x0 = new double[CW];
    double *y0 = new double[CH];
    for (i = 0; i < CW; i++)
    {
        x0[i] = x_int / 2.0 + x_int * i;
    }
    for (i = 0; i < CH; i++)
    {
        y0[i] = y_int / 2.0 + y_int * i;
    }
    double *X0 = new double[CW * CH];
    double *Y0 = new double[CW * CH];
    for (i = 0; i < CW * CH; i++)
    {
        X0[i] = x0[i - i / CW * CW];
        Y0[i] = y0[i / CW];
    }
    for (i = 0; i < CW * CH; i++)
    {
        Dxx[2 * i] = DW;
        Dxx[2 * i + 1] = DH;
        xy00[2 * i] = X0[i];
        xy00[2 * i + 1] = Y0[i];
        femxy00[2 * i] = X0[i] * DW;
        femxy00[2 * i + 1] = Y0[i] * DH;
    }
    delete x0;
    delete y0;
    delete X0;
    delete Y0;
}
