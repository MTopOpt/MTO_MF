void MaxObj_Top1(int allcells,double *hi,double *dsb,double &f,double *FF,double *Td,double *fse)
{
   
    int ii,jj,kk;
    for (ii = 0; ii < allcells; ii++)
    {
        Td[ii] = 0;
        for (jj = 0; jj < 50; jj++) 
        {
            Td[ii]=Td[ii]+dsb[jj*allcells+ii]*hi[jj];
        }
    } 

    double sumT = 0,sum=0,tmax=0,Tmax,pn=1000,T3;

    for (ii = 0; ii < allcells; ii++)
    {
        tmax = std::max(Td[ii], tmax);
    }
    MPI_Allreduce(&tmax,&Tmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    
    for (ii = 0; ii < allcells; ii++)
    {
        sum = sum + std::pow(Td[ii]/Tmax, pn);
    }

    MPI_Allreduce(&sum,&sumT,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //sumT=sumT*Tmax;

    f = Tmax*std::pow(sumT, 1.0/pn);
    T3 = std::pow(sumT, 1.0/pn-1);
    for (ii = 0; ii < 50; ii++)
    {
        fse[ii] = 0;
        for (jj = 0; jj < allcells; jj++)
        {
            fse[ii] = dsb[ii*allcells+jj] * std::pow(Td[jj]/Tmax,pn-1) + fse[ii];
        }
        fse[ii] = fse[ii] * T3;
    }
}
