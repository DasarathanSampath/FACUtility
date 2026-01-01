volVectorField U
(
    IOobject(facProps.UName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
    mesh
);

volScalarField T
(
    IOobject(facProps.TName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
    mesh
);

// --- Try to read rho (if available)
IOobject rhoHeader
(
    facProps.rhoName,
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
);

autoPtr<volScalarField> rhoPtr(nullptr);

if (rhoHeader.typeHeaderOk<volScalarField>(true))
{
    Info<< "Reading rho field: " << facProps.rhoName << nl << endl;

    rhoPtr.reset
    (
        new volScalarField
        (
            IOobject(facProps.rhoName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        )
    );
}
else
{
    if (!facProps.useThermoForRho)
    {
        FatalErrorInFunction
            << "rho not found at time " << runTime.timeName()
            << " and useThermoForRho=false. Either write rho or enable thermo fallback."
            << exit(FatalError);
    }

    Info<< "rho field not found. Creating thermo and extracting rho from thermo.rho()" << nl << endl;

    autoPtr<fluidThermo> thermo(fluidThermo::New(mesh));

    rhoPtr.reset
    (
        new volScalarField
        (
            IOobject("rho", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            thermo->rho()
        )
    );
}

const volScalarField& rho = *rhoPtr;

autoPtr<volScalarField> muEffPtr(nullptr);
autoPtr<volScalarField> nuEffPtr(nullptr);
// Detect whether muEff exists at this time
IOobject muHdr
(
    facProps.muEffName,
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
);

// Detect whether nuEff exists at this time
IOobject nuHdr
(
    "nuEff",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE
);

// A small rho floor with correct dimensions
const dimensionedScalar rhoMin("rhoMin", rho.dimensions(), SMALL);

// Prefer reading muEff if present, otherwise read nuEff if present
if (muHdr.typeHeaderOk<volScalarField>(true))
{
    muEffPtr.reset
    (
        new volScalarField
        (
            IOobject(facProps.muEffName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        )
    );
    // nuEff = muEff / rho
    nuEffPtr.reset
    (
        new volScalarField
        (
            IOobject("nuEff", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            (*muEffPtr) / max(rho, rhoMin)
        )
    );

    if (Pstream::master())
    {
        Info<< "FACUtility: Read " << facProps.muEffName
            << " and derived nuEff = muEff/rho" << nl << endl;
    }
}
else if (nuHdr.typeHeaderOk<volScalarField>(true))
{
    nuEffPtr.reset
    (
        new volScalarField
        (
            IOobject("nuEff", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        )
    );
    // muEff = nuEff * rho
    muEffPtr.reset
    (
        new volScalarField
        (
            IOobject(facProps.muEffName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
            (*nuEffPtr) * max(rho, rhoMin)
        )
    );

    if (Pstream::master())
    {
        Info<< "FACUtility: Read nuEff and derived " << facProps.muEffName
            << " = nuEff*rho" << nl << endl;
    }
}
else
{
    FatalErrorInFunction
        << "Neither " << facProps.muEffName << " nor nuEff found in time "
        << runTime.timeName() << " for region " << mesh.name() << nl
        << "Provide muEff (recommended) or nuEff in the time directory."
        << exit(FatalError);
}

// Convenience references for later use
// const volScalarField& muEff = *muEffPtr;
const volScalarField& nuEff = *nuEffPtr;

// Surface outputs (write them)
surfaceScalarField UTauFace
(
    IOobject("UTauFace", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

surfaceScalarField kmFace
(
    IOobject("kmFace", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

surfaceScalarField CeqFace
(
    IOobject("CeqFace", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

surfaceScalarField TLFace
(
    IOobject("TLFace", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

// Optional: volume-mapped for visualization
volScalarField UTau
(
    IOobject("UTau", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

volScalarField TL
(
    IOobject("TL", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

volScalarField km
(
    IOobject("km", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

volScalarField Ceq
(
    IOobject("Ceq", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

// Surface outputs as volume scalar fields only for visualization in ParaView
volScalarField UTauPatch
(
    IOobject("UTauPatch", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

volScalarField kmPatch
(
    IOobject("kmPatch", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);

volScalarField CeqPatch
(
    IOobject("CeqPatch", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

volScalarField TLPatch
(
    IOobject("TLPatch", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimVelocity, 0.0)
);