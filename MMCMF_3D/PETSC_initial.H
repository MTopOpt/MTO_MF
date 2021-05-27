//***********************************//Petsc
PetscInt n_dof=n+varnum*mylen,N_dof=N+Varnum,m=3;//m:number of constraints
PetscScalar Xmin=0; // Min. value of design variables
PetscScalar Xmax=1; // Max. value of design variables

Vec xmma,xold,xmin,xmax;
PetscScalar *xmma_array;
VecCreate(PETSC_COMM_WORLD,&xmma);
VecSetSizes(xmma,n_dof,N_dof);
VecSetFromOptions(xmma);
VecGetArray(xmma,&xmma_array);
for(i=0;i<n;i++)
{
  xmma_array[i]=x[i];
}
for (i = n; i < n_dof; i++)
{
    xmma_array[i] = xy00[i - n + mystart * varnum];
}
VecRestoreArray(xmma,&xmma_array);

VecDuplicate(xmma,&xold); 
VecDuplicate(xmma,&xmin);
VecDuplicate(xmma,&xmax); 

PetscScalar aMMA[m];
PetscScalar cMMA[m];
PetscScalar dMMA[m];

for (PetscInt i=0;i<m;i++)
{
  aMMA[i]=0.0;
  dMMA[i]=1.0;
  cMMA[i]=1000.0;
}

MMA *mma;
mma=new MMA(N_dof,m,xmma, aMMA, cMMA, dMMA,mma_init,mma_dec,mma_inc,mma_dec2,mma_inc2);

Vec dfdx; // Sensitivities of objective
VecDuplicate(xmma,&dfdx); 
PetscScalar *dfdx_array;
VecGetArray(dfdx,&dfdx_array);
VecRestoreArray(dfdx,&dfdx_array); 

PetscScalar *gx= new PetscScalar[m]; // Array with constraint values
Vec *dgdx; // Sensitivities of constraints
VecDuplicateVecs(xmma,m, &dgdx); 
PetscScalar **dgdx_array;
VecGetArrays(dgdx,m,&dgdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

