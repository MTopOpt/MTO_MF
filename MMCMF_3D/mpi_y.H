


List<scalarField> Allf2List(Pstream::nProcs());
Allf2List[Pstream::myProcNo()] = fsens2.internalField();
Pstream::gatherList(Allf2List);
Pstream::scatterList(Allf2List);
scalarField Allfsens2(ListListOps::combine<scalarField>(Allf2List,accessOp<scalarField>()));

List<scalarField> Allr2List(Pstream::nProcs());
Allr2List[Pstream::myProcNo()] = gamma2.internalField();
Pstream::gatherList(Allr2List);
Pstream::scatterList(Allr2List);
scalarField Allgamma2(ListListOps::combine<scalarField>(Allr2List,accessOp<scalarField>()));

