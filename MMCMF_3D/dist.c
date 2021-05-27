double dist(double xy[], int N)
{
    int i, j;
    double mind = 1000;
    double *dis = new double[N];
    double ln = 8 * 5e-4;

    for (i = 0; i < N; i++)
    {
        dis[i] = 1000;
        for (j = 0; j < i; j++)
        {
            dis[i] = std::min(std::max(std::abs(xy[2 * i] - xy[2 * j]), std::abs(xy[2 * i + 1] - xy[2 * j + 1])), dis[i]);
        }
        for (j = i + 1; j < N; j++)
        {
            dis[i] = std::min(std::max(std::abs(xy[2 * i] - xy[2 * j]), std::abs(xy[2 * i + 1] - xy[2 * j + 1])), dis[i]);
        }
        dis[i] = std::min(std::abs(xy[2 * i] + ln), dis[i]);
        dis[i] = std::min(std::abs(xy[2 * i] - ln - 8e-2), dis[i]);
        dis[i] = std::min(std::abs(xy[2 * i + 1] + ln), dis[i]);
        dis[i] = std::min(std::abs(xy[2 * i + 1] - ln - 0.16), dis[i]);
    }
    for (i = 0; i < N; i++)
    {
        mind = std::min(mind, dis[i]);
    }
    delete dis;
    return mind;
}
