
Info<<"sensitivity analysis"<<endl;
fsensh=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoc;
fsens2=-Tb *Q0;

gsensh_U=-alphaMax *(1 + qu) * qu / ((qu + xh) * (qu + xh)) * (U & Ua);

mu = -(1 + qd) * qd / ((qd + xh) * (qd + xh)) * (E - Emin) / (2.0 * (1.0 + Po));
lambda = -(1 + qd) * qd / ((qd + xh) * (qd + xh)) * Po * (E - Emin) / ((1.0 + Po) * (1.0 - 2.0 * Po));
gsensh_D=-gradD && (rho * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD)));

#include "filter_sens.H"

#include "MMC_sens.H"

gx[0] = PowerDiss / Foam::max(PowerDissRelax - opt * 0.2, PowerDissMax) - 1.0;
gx[1] = C / Foam::max(ComplianceRelax - opt * 0.1, ComplianceMax) - 1.0;
gx[2] = (2 * len * Ndx - mind) / len*1.1;

VecGetArray(dfdx,&dfdx_array);
VecGetArrays(dgdx,m,&dgdx_array);
for (i = 0; i < n; i++)
{
    dfdx_array[i] = fsens[i] / N;     //sensitivities of objective function
    dgdx_array[0][i] = gsens_U[i] / N; //sensitivities of constraint function
    dgdx_array[1][i] = gsens_D[i]*1.07 / N;
    dgdx_array[2][i] = 0;
}
for (i = n; i < n_dof; i++)
{
    dfdx_array[i] = Qf2[i - n]/N;
    dgdx_array[0][i] = 0;
    dgdx_array[1][i] = 0;
    dgdx_array[2][i] = Dg2[i - n]/N;
}
VecRestoreArray(dfdx,&dfdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

Info<< "\n run MMA \n" << endl;
mma->SetOuterMovelimit(Xmin,Xmax,1,movlim,xmma,xmin,xmax,n);
mma->Update(xmma,dfdx,gx,dgdx,xmin,xmax,raa0,mma_dec,mma_inc,mma_dec2,mma_inc2,n);

VecGetArray(xmma,&xmma_array);
for(i=0;i<n;i++)
{
   x[i]=xmma_array[i];
}
for (i = n; i < n_dof; i++)
{
    //femxy00[i - n + mystart * varnum] = xmma_array[i] * Dxx[i - n];
    Femxy00[i - n] = xmma_array[i] * Dxx[i - n];
}
VecRestoreArray(xmma,&xmma_array);

Allxy00[Pstream::myProcNo()] = Femxy00;
Pstream::gatherList(Allxy00);
Pstream::scatterList(Allxy00);
Allxy=ListListOps::combine<scalarField>(Allxy00,accessOp<scalarField>());
for(i=0;i<Varnum;i++)
{
    femxy00[i]=Allxy[i];
}

#include "filter_x.H"
#include "MMC_phi.H"