#include "FACErosionModel.H"
#include "fvcSnGrad.H"
#include "error.H"
#include "Pstream.H"

namespace Foam
{

FACErosionModel::FACErosionModel(const fvMesh& mesh, const FACProperties& props)
:
    mesh_(mesh), props_(props), Dtab_(nullptr)
{
    if (props_.diffusivityMode == "table")
    {
        Dtab_.reset(new FACInterpolationTable(props_.DtableFile));
    }
}

scalar FACErosionModel::diffusivity(scalar Tkelvin) const
{
    if (props_.diffusivityMode == "table" && Dtab_.valid())
        return Dtab_->value(Tkelvin);

    return props_.Dconstant;
}

/*---------------*\
 * Gradient based method
 \*--------------*/
scalar FACErosionModel::utauFromPatchLocalFace
(
    label patchId,
    label localFace,
    const surfaceVectorField& snGradU,
    const volScalarField& nuEff,
    const volVectorField& U
) const
{
    const fvPatch& p = mesh_.boundary()[patchId];
    const vectorField nHat(p.nf());
    const vector n = nHat[localFace];

    const vector dUdn = snGradU.boundaryField()[patchId][localFace];

    const vector dUdnTan = dUdn - (dUdn & n)*n;     
    const scalar magGradTan = mag(dUdnTan);

    const scalar nuE = nuEff.boundaryField()[patchId][localFace];
    const scalar UtauGradientMethod = sqrt(nuE * magGradTan);
    return UtauGradientMethod;
}

void FACErosionModel::compute
(
    const FACSelection& sel,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& nuEff,
    surfaceScalarField& UTauFace,
    surfaceScalarField& kmFace
) const
{
    // Initialize to zero (for non-selected faces)
    UTauFace = dimensionedScalar("zero", dimVelocity, 0.0);
    kmFace   = dimensionedScalar("zero", dimVelocity, 0.0);

    const labelList& faces = sel.globalFaces();

    // snGrad(U) gives dU/dn vector on faces
    const surfaceVectorField snGradU = fvc::snGrad(U);

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

        // --- 1) utau on that boundary face
        const scalar utau = utauFromPatchLocalFace(pId, lf, snGradU, nuEff, U);

        // store in UTauFace (boundary face)
        UTauFace.boundaryFieldRef()[pId][lf] = utau;

        const label owner = mesh_.faceOwner()[gf];
        const scalar Um = props_.UmeanConstant;

        const scalar nuEcell = nuEff.boundaryField()[pId][lf];
        const scalar Re = Um * props_.charLength / nuEcell;

        const scalar D = diffusivity(T[owner]);
        const scalar Sc = nuEcell / max(D, SMALL);

        // friction factor: f/2 = (utau/Um)^2 (Eq.15 equivalent)
        const scalar fOver2 = sqr(utau / max(Um, SMALL));

        const scalar km = fOver2 * Re * pow(Sc, 0.33) * D / max(props_.charLength, SMALL);

        kmFace.boundaryFieldRef()[pId][lf] = km;
    }
}

} // namespace Foam
