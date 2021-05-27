if (opt % 10 == 0||opt==1) //
{
    if (Pstream::master())
{ 
   ofstream outfile24("Time2.txt",std::ios::app);
   outfile24<<"exchange_begain="<<runTime.elapsedClockTime()<<"\n";
   outfile24.close();   
}
    for (i = 0; i < 50; i++)
    {
        HI[i] = hi[i];
    }
    for (ii = 0; ii < 50; ii++)
    {
        for (jj = 0; jj < 50; jj++)
        {
            hi[jj]=0;
        }
        hi[ii] = 1;
        #include "MMC_phi.H"
        #include "Primal_T.H"
        for (k = 0; k < n; k++)
        {
            dsb[ii * n + k] = T[k];
        }
    }
    if (Pstream::master())
{ 
   ofstream outfile24("Time2.txt",std::ios::app);
   outfile24<<"T_end="<<runTime.elapsedClockTime()<<"\n";
   outfile24.close();   
}
    while (1)
    {
        MaxObj_Top1(n, HI, dsb, f, Td, fse);//f:maxT,fse:sensitivity,Td:T
 
        if (Pstream::master())
        {
            ofstream outfile444("ffff.txt", std::ios::app);
            outfile444 << "opt=" << opt << "f0=" << f << "\n";
            outfile444.close();
        }
        MPI_Allreduce(&fse, &FF, 50, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (i = 0; i < 50; i++)
        {
            for (j = 0; j < 50; j++)
            {
                KFQ[i * 50 + j] = (HI[i] - HI[j]) * (FF[i] - FF[j]);
            }
        }

        while (1)
        {
            for (k = 0; k < 50; k++)
            {
                hhi[k] = HI[k];
            }
            sort(KFQ, imax, jmax, KFQmax);
            if (KFQmax <= 0)
            {
                break;
            }
            tem = hhi[imax];
            hhi[imax] = hhi[jmax];
            hhi[jmax] = tem;
            KFQ[imax * 50 + jmax] = 0;
            KFQ[jmax * 50 + imax] = 0;
            MaxObj_Top2(allcells, hhi, dsb, f1, Td);
            if (f1 < f)
            {
                for (k = 0; k < 50; k++)
                {
                    HI[k] = hhi[k];
                }
                if (Pstream::master())
                {
                    ofstream outfile444("ffff.txt", std::ios::app);
                    outfile444 << "opt=" << opt << "f1=" << f1 << "\n";
                    outfile444.close();
                }
                break;
            }
        }
        sort(KFQ, imax, jmax, KFQmax);
        if (KFQmax <= 0)
        {
            break;
        }
    }
    if (Pstream::master())
{ 
   ofstream outfile24("Time2.txt",std::ios::app);
   outfile24<<"exchange_end="<<runTime.elapsedClockTime()<<"\n";
   outfile24.close();   
}
    for (i = 0; i < 50; i++)
    {
        hi[i] = HI[i];
    }
    #include "MMC_phi.H"
}