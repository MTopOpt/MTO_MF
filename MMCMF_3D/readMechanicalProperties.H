Info<< "Reading mechanical properties\n" << endl;
IOdictionary mechanicalProperties
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
int iCorr=0,nCorr=300;
scalar convergenceTolerance(5e-5),initialResidual(0);

Switch planeStress(mechanicalProperties.lookup("planeStress"));
dimensionedScalar rhoE("rhoE", dimMass/dimLength/sqr(dimTime),readScalar(mechanicalProperties.lookup("rhoE")));
dimensionedScalar rho("rho",  dimMass/dimVolume,readScalar(mechanicalProperties.lookup("rho")));
scalar Po(readScalar(mechanicalProperties.lookup("Po")));
dimensionedScalar E(rhoE/rho);
dimensionedScalar Emin(E*1e-4);
volScalarField mu(xh*xh*xh*E/(2.0*(1.0 + Po)));
volScalarField lambda(xh*xh*xh*Po*E/((1.0 + Po)*(1.0 - 2.0*Po)));
scalar qd(readScalar(mechanicalProperties.lookup("qd")));
Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(D)
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(gradD) + (lambda*I)*tr(gradD)
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

Switch compactNormalStress("yes");

if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());

volScalarField gsens_D
(
    IOobject
    (
        "gsens_D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -gradD && (rho * (mu * twoSymm(gradD) + (lambda * I) * tr(gradD))),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsensh_D(gsens_D);

volSymmTensorField sigma(rho *(mu *twoSymm(gradD) + (lambda * I) * tr(gradD)));
volScalarField Energy((gradD + gradD.T()) && sigma);