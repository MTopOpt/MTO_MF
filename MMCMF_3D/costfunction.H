Info<<"cost function"<<endl;

gradD = fvc::grad(D);
sigma=rho *(mu *twoSymm(gradD) + (lambda * I) * tr(gradD));
Energy=(gradD + gradD.T()) && sigma;
C=fvc::domainIntegrate(Energy).value()/Compliance0/2.0;
V=fvc::domainIntegrate(xh).value()/area;  
sumQ = Quse - gSum(gamma2)/Q_dof;

PowerDiss=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDiss=PowerDiss/PowerDiss0;
if(GeoDim==2)
{
 PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>()); 

if (Pstream::master())
{ 
   ofstream outfile1("Voluse.txt",std::ios::app);
   outfile1<<V<<"\n";
   outfile1.close();

   ofstream outfile2("PowerDiss.txt",std::ios::app);
   outfile2<<PowerDiss<<"\n";
   outfile2.close(); 

   ofstream outfile3("sumQ.txt",std::ios::app);
   outfile3<<Quse -sumQ<<"\n";
   outfile3.close(); 

   ofstream outfile4("Time.txt",std::ios::app);
   outfile4<<runTime.elapsedClockTime()<<"\n";
   outfile4.close();    

   ofstream outfile5("maxT.txt",std::ios::app);
   outfile5<<maxT<<"\n";
   outfile5.close();

    ofstream outfile6("compliance.txt", std::ios::app);
    outfile6 << C << "\n";
    outfile6.close();

    ofstream outfile33("femxy00.txt", std::ios::app);
    outfile33 << opt << "\n";
    for (i = 0; i < Varnum; i++)
    {
        outfile33 << femxy00[i] << "\n";
    }
    outfile33 << "\n";
    outfile33.close();
}
if (runTime.writeTime())
{
    T.write();
}