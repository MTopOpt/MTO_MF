Info << "MMC_phi" << endl;

//solve TDF of each components
for (i = mystart; i <=myend; i++)
{
    for (j = 0; j < varnum; j++)
    {
        Cxy00[j] = femxy00[varnum * i + j];
    }

    MMC_tPhi(Cxy00, dx, dy, cellpos_x, cellpos_z, Q_dof, tphi, pTDF, hi[i]);

    for (j = 0; j < Q_dof; j++)
    {
        PHI[(i - mystart) * Q_dof + j] = tphi[j];
    }
}

//
AllPHI[Pstream::myProcNo()] = PHI;
Pstream::gatherList(AllPHI);
Pstream::scatterList(AllPHI);
Allphi=ListListOps::combine<scalarField>(AllPHI,accessOp<scalarField>());

//
for (j = 0; j < Q_dof; j++)
{
    Allgamma2[Q_index[j]] = Allphi[j];
}

for (i = 1; i < Comnum; i++)
{
    for (j = 0; j < Q_dof; j++)
    {
        Allgamma2[Q_index[j]] = Foam::max(Allphi[i * Q_dof + j], Allgamma2[Q_index[j]]);
    }
}

for (i = 0; i < n; i++)
{
    gamma2[i] = Allgamma2[i + Ce];
}
gamma2.correctBoundaryConditions();
