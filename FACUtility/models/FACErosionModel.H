#pragma once
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "FACSelection.H"
#include "FACInterpolation.H"
#include "FACProperties.H"

namespace Foam
{

class FACErosionModel
{
    const fvMesh& mesh_;
    const FACProperties& props_;
    autoPtr<FACInterpolationTable> Dtab_;

public:
    FACErosionModel(const fvMesh& mesh, const FACProperties& props);

    void compute
    (
        const FACSelection& sel,
        const volVectorField& U,
        const volScalarField& T,
        const volScalarField& nuEff,
        surfaceScalarField& UTauFace,
        surfaceScalarField& kmFace
    ) const;

private:
    scalar diffusivity(scalar Tkelvin) const;

    scalar utauFromPatchLocalFace
    (
        label patchId,
        label localFace,
        const surfaceVectorField& snGradU,
        const volScalarField& nuEff,
        const volVectorField& U
    ) const;
};

}
