Info<< "Reading thermal properties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoc(thermalProperties.lookup("rhoc"));//rho_fluid*heat_capacity
dimensionedScalar kf(thermalProperties.lookup("kf"));//Thermal Conductivity of fluid
dimensionedScalar ks(thermalProperties.lookup("ks"));//Thermal Conductivity of solid
//dimensionedScalar Q("Q", dimensionSet(0,0,-1,1,0,0,0),1.0);
dimensionedScalar b2("b2", dimensionSet(0,0,-1,2-Pnorm,0,0,0),1);
  volScalarField DT
   (
       IOobject
       (
           "DT",
           runTime.timeName(),
           mesh,
           IOobject::READ_IF_PRESENT,
           IOobject::AUTO_WRITE
       ),
       (kf+(ks-kf)*x*x*x)/rhoc
   );

Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
Info<< "Reading field Tb\n" << endl;
   volScalarField Tb
   (
       IOobject
       (
           "Tb",
           runTime.timeName(),
           mesh,
           IOobject::MUST_READ,
           IOobject::AUTO_WRITE
       ),
       mesh
   ); 
volScalarField T2(T);

     volScalarField Q0
     (
         IOobject
         (
             "Q0",
             runTime.timeName(),
             mesh,
             IOobject::READ_IF_PRESENT,
             IOobject::NO_WRITE
         ),
         mesh,
         dimensionedScalar("Q0", dimensionSet(0,0,-1,1,0,0,0),200000000.0/4180000),
         calculatedFvPatchScalarField::typeName
     );
     
volScalarField fsens
(
    IOobject
    (
        "fsens",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField fsensh(fsens);
volScalarField fsens2(-Tb *Q0);