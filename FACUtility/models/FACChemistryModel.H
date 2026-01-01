#pragma once

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "autoPtr.H"

#include "FACInterpolation.H"
#include "FACChemistryConstants.H"
#include "FACNewtonRaphson.H"
#include "FACSelection.H"
#include "FACProperties.H"   // <-- your FACProperties struct/reader

namespace Foam
{

class FACChemistryModel
{
    const fvMesh& mesh_;
    const FACProperties& props_;

    autoPtr<FACInterpolationTable> pH2Tab_;
    autoPtr<FACInterpolationTable> KwTab_;
    ChemistryConstants chem_;

    scalar mA_;                 // KOH molality (from FACProperties)
    NewtonControls newton_;     // can be filled from FACProperties too

    // Kb(T) for each species
    scalar Kb(const SpeciesConst& sp, const scalar T) const;

    scalar FmH
    (
        const scalar mH,
        const scalar k0, const scalar k1, const scalar k2, const scalar k3,
        const scalar pH2, const scalar Kw, const scalar mA
    ) const;

    scalar dFmH
    (
        const scalar mH,
        const scalar k0, const scalar k1, const scalar k2, const scalar k3,
        const scalar pH2, const scalar Kw, const scalar mA
    ) const;

public:
    FACChemistryModel(const fvMesh& mesh, const FACProperties& facProps);

    void compute
    (
        const FACSelection& sel,
        const volScalarField& T,
        surfaceScalarField& CeqFace
    ) const;
};

} // namespace Foam
