
VecDestroy(&xmma);
VecDestroy(&xold);
VecDestroy(&xmin);
VecDestroy(&xmax);
VecDestroy(&dfdx);
VecDestroyVecs(m, &dgdx);
delete[] gx;
delete mma;

delete Dxx ; 
delete femxy00 ; 
delete Cxy00 ; 
delete tphi ; 
delete tphi2 ; 
delete tphi3 ; 
delete cellpos_x ; 
delete cellpos_z ; 
delete Q_index ; 
delete q_index ; 

Info<< " \n" << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
<< nl << endl;