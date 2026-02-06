// -------------------------------------------------------------------------
// Primary flow fields
// -------------------------------------------------------------------------

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field U_old\n" << endl;
volVectorField U_old
(
    IOobject
    (
        "U_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U
);

// -------------------------------------------------------------------------
// MHD coupling fields
// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// MHD coupling fields (Elmer -> OpenFOAM)
// Elmer exports Volume Current as components: Jx, Jy, Jz
// We reconstruct J as a volVectorField
// -------------------------------------------------------------------------

Info<< "Reading magnetic field B\n" << endl;

volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing current density components (Jx, Jy, Jz)\n" << endl;

volScalarField Jx
(
    IOobject
    (
        "Jx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent/dimArea,   
        0.0
    )
);

volScalarField Jy
(
    IOobject
    (
        "Jy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent/dimArea,
        0.0
    )
);

volScalarField Jz
(
    IOobject
    (
        "Jz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimCurrent/dimArea,
        0.0
    )
);

Info<< "Initializing reconstructed current vector J (from components)\n" << endl;

volVectorField J_dens
(
    IOobject
    (
        "J",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector
    (
        "zero",
        dimCurrent/dimArea,
        vector::zero
    )
);

Info<< "Initializing Lorentz force field fLorentz\n" << endl;

volVectorField fLorentz
(
    IOobject
    (
        "fLorentz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    J_dens ^ B // This is cross product for openFOAM
);



Info<< "Initializing field elcond\n" << endl;
volScalarField elcond
(
    IOobject
    (
        "elcond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimensionSet(-1, -3, 3, 0, 0, 2, 0),
        0.0
    )
);

Info<< "Initializing field T\n" << endl;
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

Info<< "Initializing field JH\n" << endl;
volScalarField JH
(
    IOobject
    (
        "JH",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Initializing field JH_recv\n" << endl;
volScalarField JH_recv
(
    IOobject
    (
        "JH_recv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    JH
);

// -------------------------------------------------------------------------
// Fluxes and transport
// -------------------------------------------------------------------------

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
singlePhaseTransportModel laminarTransport(U, phi);

// Constant fluid properties
dimensionedScalar rhoConst("rho", dimDensity, laminarTransport);
dimensionedScalar cp("cp", dimSpecificHeatCapacity, laminarTransport);
dimensionedScalar kappa("kappa", dimPower/dimLength/dimTemperature, laminarTransport);

// Density field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    rhoConst
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

// Turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// -------------------------------------------------------------------------
// Gravity and pressure
// -------------------------------------------------------------------------

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());

// -------------------------------------------------------------------------
// MRF (after mesh & transport are ready)
// -------------------------------------------------------------------------

#include "createMRF.H"

// -------------------------------------------------------------------------
// Electrical conductivity (single-body enforced)
// -------------------------------------------------------------------------

dimensionedScalar elcond_melt
(
    "elcond_melt",
    dimensionSet(-1, -3, 3, 0, 0, 2, 0),
    laminarTransport
);

// Force single, nonzero conductivity everywhere
elcond = elcond_melt;

// -------------------------------------------------------------------------
// Elmer update control
// -------------------------------------------------------------------------

// -------------------------------------------------------------------------
// Scalar velocity components for Elmer coupling
// -------------------------------------------------------------------------

Info<< "Initializing velocity components (Ux, Uy, Uz)\n" << endl;

volScalarField Ux
(
    IOobject
    (
        "Ux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimVelocity,
        0.0
    )
);

volScalarField Uy
(
    IOobject
    (
        "Uy",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimVelocity,
        0.0
    )
);

volScalarField Uz
(
    IOobject
    (
        "Uz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimVelocity,
        0.0
    )
);

// Initialize once
Ux = U.component(vector::X);
Uy = U.component(vector::Y);
Uz = U.component(vector::Z);


// -------------------------------------------------------------------------
// Scalar magnetic field components for Elmer coupling
// -------------------------------------------------------------------------

Info<< "Initializing magnetic field components (Bx, By, Bz)\n" << endl;

volScalarField Bx
(
    IOobject
    (
        "Bx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        B.dimensions(),
        0.0
    )
);

volScalarField By
(
    IOobject
    (
        "By",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        B.dimensions(),
        0.0
    )
);

volScalarField Bz
(
    IOobject
    (
        "Bz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        B.dimensions(),
        0.0
    )
);

// Initialize once
Bx = B.component(vector::X);
By = B.component(vector::Y);
Bz = B.component(vector::Z);

