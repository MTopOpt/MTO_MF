//***********************************//update mu and lambda
opt++;

if(solid_area)
{
   forAll(cells_solid, i)
   {
      xh[cells_solid[i]]=0;
      x[cells_solid[i]]=0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      xh[cells_fluid[i]]=1.0;
      x[cells_fluid[i]]=1.0;
   }
}
if(Q_area)
{
   forAll(cells_Q, i)
   {
      xh[cells_Q[i]]=0;
      x[cells_Q[i]]=0;

   }
}
xh.correctBoundaryConditions();
if(runTime.writeTime())
{
   xh.write();  
   gamma2.write();
}


if(opt>=100)//45,95,145
{
  alphaMax=alphaMax+alphamax/30.0;
}

if(opt>=80)
{
  qu=0.005+(opt-80)*1e-4;
}
qu = Foam::min(qu, 0.01);
alpha = alphaMax * qu * (1 - xh) / (qu + xh);
DT = (ks + (kf - ks) * xh * (1 + qu) / (qu + xh)) / rhoc;

mu = qd * (1 - xh) / (qd + xh) * (E - Emin) / (2.0 * (1.0 + Po)) + Emin / (2.0 * (1.0 + Po));
lambda = qd * (1 - xh) / (qd + xh) * Po * (E - Emin) / ((1.0 + Po) * (1.0 - 2.0 * Po)) + Po * Emin / ((1.0 + Po) * (1.0 - 2.0 * Po));

mind = dist(femxy00, Comnum);
if (Pstream::master())
{
    ofstream outfile12("mind.txt", std::ios::app);
    outfile12 << "iter=" << opt << "\n";
    outfile12 << "mind=" << mind / len << "\n";
    outfile12 << "\n";
    outfile12.close();
}