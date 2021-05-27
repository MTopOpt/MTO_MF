void sort(double *KFQ,int &imax,int &jmax,double &KFQmax)
{
    int i,j;
    KFQmax=-10000;
    for(i=0;i<50;i++)
    {
        for(j=0;j<50;j++)
        {
          if(KFQmax<KFQ[i*50+j])
          {
              imax=i;
              jmax=j;
              KFQmax=KFQ[i*50+j];
          }
        }
    }
}