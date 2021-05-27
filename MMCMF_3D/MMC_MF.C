//Author: Yu Minghao    Updated: May 2020 
static char help[] = "topology optimization of thermal-fluid problem\n";
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"//
#include <mpi.h>
#include "Head.h"//
#include <MMA.h>

void MaxObj_Top1(int allcells,double *hi,double *dsb,double &f,double *Td,double *fse);
void MaxObj_Top2(int allcells,double *hi,double *dsb,double &f,double *Td);
int main(int argc, char *argv[])
{
    PetscInitialize(&argc,&argv,PETSC_NULL,help);
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"//
    #include "createFields.H"
    #include "readTransportProperties.H" 
    #include "initContinuityErrs.H"
    #include "readMechanicalProperties.H"
    #include "readThermalProperties.H" 
    #include "mpi_y.H"
    #include "opt_initialization.H"
    while (simple.loop(runTime))
    {
        #include "update.H"
        #include "Primal_U.H"
        #include "Exchange.H"
        #include "Primal_T.H"
        #include "Primal_D.H"
        #include "AdjointHeat_Tb.H"
        #include "AdjointHeat_Ub.H"
        #include "AdjointFlow_Ua.H"
        #include "costfunction.H"              
        #include "sensitivity.H"
    }
    #include "finalize.H"
    return 0;
}
void MaxObj_Top1(int allcells,double *hi,double *dsb,double &f,double *Td,double *fse)
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
void MaxObj_Top2(int allcells,double *hi,double *dsb,double &f,double *Td)
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
}


