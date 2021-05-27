Info << "\nMMC_sens\n" << endl;

Allf2List[Pstream::myProcNo()] = fsens2.internalField();
Pstream::gatherList(Allf2List);
Pstream::scatterList(Allf2List);
Allfsens2=ListListOps::combine<scalarField>(Allf2List,accessOp<scalarField>());

for (i = mystart; i <=myend; i++)
{
    for (j = 0; j < varnum; j++)
    {
        Cxy00[j] = femxy00[varnum * i + j];
    }
    for (j = 0; j < varnum; j++)
    {
        Cxy00[j] = Cxy00[j] + len;
        MMC_tPhi(Cxy00, dx, dy, cellpos_x, cellpos_z, Q_dof, tphi2, pTDF, hi[i]);
        Cxy00[j] = Cxy00[j] - len;

        for (ik = 0; ik < i; ik++)
        {
            for (jk = 0; jk < Q_dof; jk++)
            {
                tphi2[jk] = Foam::max(Allphi[ik * Q_dof + jk], tphi2[jk]);
            }
        }
        for (ik = i + 1; ik < Comnum; ik++)
        {
            for (jk = 0; jk < Q_dof; jk++)
            {
                tphi2[jk] = Foam::max(Allphi[ik * Q_dof + jk], tphi2[jk]);
            }
        }
        femxy00[i*varnum+j]=femxy00[i*varnum+j]+len;
        mind2=dist(femxy00, Comnum);
        femxy00[i*varnum+j]=femxy00[i*varnum+j]-len;
        //************************************//
        Cxy00[j] = Cxy00[j] - len;
        MMC_tPhi(Cxy00, dx, dy, cellpos_x, cellpos_z, Q_dof, tphi3, pTDF, hi[i]);
        Cxy00[j] = Cxy00[j] + len;
 
        for (ik = 0; ik < i; ik++)
        {
            for (jk = 0; jk < Q_dof; jk++)
            {
                tphi3[jk] = Foam::max(Allphi[ik * Q_dof + jk], tphi3[jk]);
            }
        }
        for (ik = i + 1; ik < Comnum; ik++)
        {
            for (jk = 0; jk < Q_dof; jk++)
            {
                tphi3[jk] = Foam::max(Allphi[ik * Q_dof + jk], tphi3[jk]);
            }
        }
        femxy00[i*varnum+j]=femxy00[i*varnum+j]-len;
        mind3=dist(femxy00, Comnum);
        femxy00[i*varnum+j]=femxy00[i*varnum+j]+len;
        //************************************//
        Qf2[(i - mystart) * varnum + j] = 0;
        Qg2[(i - mystart) * varnum + j] = 0;
        for (k = 0; k < Q_dof; k++)
        {
            Qf2[(i - mystart) * varnum + j] = Qf2[(i - mystart) * varnum + j] + Allfsens2[Q_index[k]] * (tphi2[k] - tphi3[k]) / (2 * delta[j])*4000;
        }
        Dg2[(i - mystart) * varnum + j]=-(mind2-mind3)/ (2 * delta[j]);
    } //j
} //i
if (gMax(mag(Dg2)) > 0.000001)
{
    Dg2 = Dg2 / gMax(mag(Dg2)) * gMax(mag(Qf2)) * 100;
}
/*
ofstream outfile132("fsens.txt", std::ios::app);
outfile132 << "iter=" << opt << "\n";
for (i = 0; i < mylen*varnum; i++)
{
    outfile132 << "qf2=" << Qf2[i] << "  "<<"dg2=" << Dg2[i] << "\n";
}
outfile132.close();
*/
