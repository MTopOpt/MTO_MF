
Maxf=gMax(mag(fsensh.primitiveField()));
fsensh=fsensh/Maxf;
fsens2=fsens2/Maxf;
gsensh_U=gsensh_U/gMax(mag(gsensh_U.primitiveField()));
gsensh_D=gsensh_D/gMax(mag(gsensh_D.primitiveField()));
//gsenshVol.primitiveFieldRef()=mesh.V()/gMax(mesh.V());

//***********************************//Heaviside_filter
if(opt>1)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-xp[i]/eta5))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+Foam::exp(-del); 
     } 
  }
  fsensh=fsensh*drho;
  gsensh_U=gsensh_U*drho;
  gsensh_D=gsensh_D*drho;
  //gsenshVol=gsenshVol*drho;
}
//***********************************//PDE filter

solve(fvm::laplacian(fsens)-fvm::Sp(b,fsens)+fsensh*b);

solve(fvm::laplacian(gsens_U)-fvm::Sp(b,gsens_U)+gsensh_U*b);

solve(fvm::laplacian(gsens_D)-fvm::Sp(b,gsens_D)+gsensh_D*b);

Maxf=gMax(mag(fsens.primitiveField()));
fsens=fsens/Maxf;
fsens2=fsens2/Maxf;
gsens_U=gsens_U/gMax(mag(gsens_U.primitiveField()));
gsens_D=gsens_D/gMax(mag(gsens_D.primitiveField()));