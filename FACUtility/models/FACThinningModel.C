#include "FACThinningModel.H"
#include "error.H"

namespace Foam
{

FACThinningModel::FACThinningModel(const fvMesh& mesh, const FACProperties& props)
:
    mesh_(mesh), 
    props_(props),
    chromiumFactor_(props_.chromiumFactor),
    corrosionHours_(props_.corrosionHours),
    chem_()
{
    chem_.read(props_.chemistryTableFile);
}

void FACThinningModel::compute
(
    const FACSelection& sel,
    const volScalarField& rho,
    const surfaceScalarField& kmFace,
    const surfaceScalarField& CeqFace,
    surfaceScalarField& TLFace
) const
{
    Info << "\n Entered TL Compute\n" << endl;

    TLFace = dimensionedScalar("zero", dimVelocity, 0.0);

    const labelList& faces = sel.globalFaces();

    const scalar Molor_mass_iron = chem_.MFe;
    const scalar Density_iron = chem_.rhoFe;

    forAll(faces, i)
    {
        const label gf = faces[i];

        if (!sel.isBoundaryFace(i))
        {
            if (props_.skipInternalFaces) continue;
            FatalErrorInFunction << "Internal face in selection. FAC requires wall boundary faces."
                                 << exit(FatalError);
        }

        const label pId = sel.patchId(i);
        const label lf  = sel.patchLocalFace(i);

        const label owner = mesh_.faceOwner()[gf];
        const scalar rhoC = mag(rho[owner]);

        const scalar km  = kmFace.boundaryField()[pId][lf];
        const scalar Ceq = CeqFace.boundaryField()[pId][lf];
        
        //Find ferrous iron solubility Ceq in L/L
        const scalar CeqLbyL = Ceq * (Molor_mass_iron/Density_iron) * rhoC * 1e-3;

        TLFace.boundaryFieldRef()[pId][lf] = chromiumFactor_ * corrosionHours_ * km * CeqLbyL * 60*60*1000;

        if (i == 0 && Pstream::master())
        {
            Info<< "FAC Thinning constants (cell " << owner << ")" << nl
                <<"rho=" << rhoC << endl;
        }
    }
    Info << "\n Exit TL Compute\n" << endl;
}

} // namespace Foam
