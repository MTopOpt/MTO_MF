Info<<"MMC_initial"<<endl;

int CW=readScalar(optProperties.lookup("CW"));
int CH=readScalar(optProperties.lookup("CH"));
double Ndx=readScalar(optProperties.lookup("Ndx"));
double Ndz=readScalar(optProperties.lookup("Ndz"));
double pTDF=readScalar(optProperties.lookup("pTDF"));

int Comnum=CW*CH;
int varnum=2,Varnum=varnum*Comnum;
int mycom=Pstream::myProcNo()*Comnum/Pstream::nProcs();
int mystart,myend,mylen;
if(Pstream::myProcNo()<10)
{
    mystart=2*Pstream::myProcNo();
    myend=mystart+1;
    mylen=2;
}
else
{
    mystart=20+3*(Pstream::myProcNo()-10);
    myend=mystart+2;
    mylen=3;
}

double DW=8e-2;//
double DH=1.6e-1;//
double dx=Ndx*len;
double dy=Ndz*len;
//double pp=5;
double mind=0,mind2,mind3;
//
double *Dxx =new double[Varnum]; 
double *xy00 =new double[Varnum]; 
double *femxy00 =new double[Varnum]; 
double *Cxy00 =new double[varnum]; 

double *tphi =new double[Q_dof]; 
double *tphi2=new double[Q_dof];
double *tphi3=new double[Q_dof];

double delta[2];
delta[0]=len/DW;
delta[1]=len/DH;

Info<<"MMC_initial"<<endl;
MMC_initial(CW,CH,DW,DH,Dxx,xy00,femxy00);

double hi[50];
double data;
ifstream dfile("Q.txt");
for(i=0;i<50;i++)
{
    dfile >> data;
	hi[i]=data;
}
dfile.close();
int imax,jmax,ii,jj;
double KFQ[2500];
double f=0,f1,KFQmax,tem;
double FF[50],hhi[50],fse[50],HI[50];

double *Td = new double [n];
double *dsb = new double [50*n];
for(i=0;i<50;i++)
{
    FF[i]=0;
}

scalarField PHI(mylen*Q_dof,-1);
scalarField Qf2(mylen*varnum,0);
scalarField Qg2(mylen*varnum,0);
scalarField Dg2(mylen*varnum,0);

List<scalarField> AllPHI(Pstream::nProcs());
AllPHI[Pstream::myProcNo()] = PHI;
Pstream::gatherList(AllPHI);
Pstream::scatterList(AllPHI);
scalarField Allphi(ListListOps::combine<scalarField>(AllPHI,accessOp<scalarField>()));

scalarField Femxy00(varnum*mylen,-1);
List<scalarField> Allxy00(Pstream::nProcs());
Allxy00[Pstream::myProcNo()] = Femxy00;
Pstream::gatherList(Allxy00);
Pstream::scatterList(Allxy00);
scalarField Allxy(ListListOps::combine<scalarField>(Allxy00,accessOp<scalarField>()));


#include "MMC_phi.H"





