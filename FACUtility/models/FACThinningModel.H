#pragma once
#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "FACSelection.H"
#include "FACProperties.H"

#include "FACChemistryConstants.H"
#include "autoPtr.H"

namespace Foam
{

class FACThinningModel
{
    const fvMesh& mesh_;
    const FACProperties& props_;

    scalar chromiumFactor_;
    scalar corrosionHours_;

    ChemistryConstants chem_;

public:
    FACThinningModel(const fvMesh& mesh, const FACProperties& props);

    void compute
    (
        const FACSelection& sel,
        const volScalarField& rho,
        const surfaceScalarField& kmFace,
        const surfaceScalarField& CeqFace,
        surfaceScalarField& TLFace
    ) const;
};

}
