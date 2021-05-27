Info<< "Initialization\n" << endl;
std::remove("Voluse.txt");
std::remove("Time.txt");
std::remove("sumQ.txt"); 
std::remove("maxT.txt");
std::remove("PowerDiss.txt"); 
std::remove("compliance.txt"); 
std::remove("fsens.txt"); 
std::remove("femxy00.txt"); 
std::remove("mind.txt"); 
std::remove("Time2.txt"); 

int i=0,j=0,k=0,ik,jk,opt=0,Ce=0;
double del=0;
scalar allcells(0),nallcells(0),eta0,eta5,eta1,y0,y5,Maxf(0);
allcells=x.size();//n:length of x in the i-th core
nallcells=x.size();//N:design freedom
reduce(nallcells, sumOp<scalar>()); 
int n=allcells,N=nallcells;

List<scalar> cellsize(Pstream::nProcs());
for (i = 0; i < Pstream::nProcs(); i++)
{
    cellsize[i] = 0;
}
cellsize[Pstream::myProcNo()] = n;
reduce(cellsize, sumOp<List<scalar>>());
for (i = 0; i < Pstream::myProcNo(); i++)
{
    Ce = Ce + cellsize[i];
}

scalar C(0),V(0),area=gSum(mesh.V()),PowerDiss(0),maxT(0),sumT(0),sumQ(0);
double len=5e-4;
b.value()=1.0/(filterR*len/3.464)/(filterR*len/3.464);
//mesh
List<vectorField> AllcellList(Pstream::nProcs());
vectorField localcell(mesh.C());
AllcellList[Pstream::myProcNo()] = localcell;
Pstream::gatherList(AllcellList);
Pstream::scatterList(AllcellList);
vectorField Allcellcenter(
    ListListOps::combine<vectorField>(
        AllcellList,
        accessOp<vectorField>()));

int Q_dof = 0, q_dof = 0;
for (i = 0; i < N; i++)
{
    if (Allcellcenter[i][1] > (8e-3 - len))
    {
        Q_dof++;
    }
}
for (i = 0; i < n; i++)
{
    if (mesh.C()[i][1] > (8e-3 - len))
    {
        q_dof++;
    }
}
double *cellpos_x = new double[Q_dof];
double *cellpos_z = new double[Q_dof];
int *Q_index = new int[Q_dof];
int *q_index = new int[q_dof];
j = 0;
for (i = 0; i < n; i++)
{
    if (mesh.C()[i][1] > (8e-3 - len))
    {
        q_index[j] = i;
        j++;
    }
}
j = 0;
for (i = 0; i < N; i++)
{
    if (Allcellcenter[i][1] > (8e-3 - len))
    {
        cellpos_x[j] = Allcellcenter[i][0];
        cellpos_z[j] = Allcellcenter[i][2];
        Q_index[j] = i;
        j++;
    }
}
#include "MMC_initial.H"
#include "PETSC_initial.H"