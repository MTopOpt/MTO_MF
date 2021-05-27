double dist(double xy[], int N)
{
    int i, j;
    double *dis = new double[N], mind = 0;
    for (i = 0; i < N; i++)
    {
        dis[i] = 1000;
        for (j = 0; j < i; j++)
        {
            dis[i] = std::min(std::max(std::abs(xy[i] - xy[j]), std::abs(xy[i + 1] - xy[j + 1])), dis[i]);
        }
        for (j = i + 1; j < N; j++)
        {
            dis[i] = std::min(std::max(std::abs(xy[i] - xy[j]), std::abs(xy[i + 1] - xy[j + 1])), dis[i]);
        }
    }
    for (i = 0; i < N; i++)
    {
        mind=std::min(mind,dis[i]);
    }
    return mind;
}
