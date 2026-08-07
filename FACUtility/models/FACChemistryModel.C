#include "FACChemistryModel.H"
#include "error.H"
#include "mathematicalConstants.H"
#include "Pstream.H"

namespace Foam
{

FACChemistryModel::FACChemistryModel
(
    const fvMesh& mesh,
    const FACProperties& props
)
:
    mesh_(mesh),
    props_(props),
    pH2Tab_(new FACInterpolationTable(props_.pH2TableFile)),
    KwTab_(new FACInterpolationTable(props_.KwTableFile)),
    chem_(), // or chemistryConstantsFile
    mA_(props_.mA),
    newton_(props_.newton)
{
    chem_.read(props_.chemistryTableFile);

    if (chem_.species.size() < 4)
    {
        FatalErrorInFunction
            << "chemistryConstants file must provide at least 4 speciesConstants rows. "
            << "Found: " << chem_.species.size()
            << exit(FatalError);
    }

    // Print what got loaded once
    if (Pstream::master())
    {
        Info<< "FAC chemistry: Rcal=" << chem_.Rcal
            << ", species rows=" << chem_.species.size()
            << ", pH2Table=" << props_.pH2TableFile
            << ", KwTable=" << props_.KwTableFile
            << ", chemFile=" << props_.chemistryTableFile
            << nl << endl;
    }

    if (props_.pH2TableFile.empty() || props_.KwTableFile.empty() || props_.chemistryTableFile.empty())
    {
        FatalErrorInFunction << "Missing chemistry files in system/FACProperties" << exit(FatalError);
    }

}
// ---- IMPORTANT: Put your VERIFIED Kb(T) equation here ----
// You have Ab [cal/mol], Bb [cal/mol-K], Db [cal/mol-K], R [cal/mol-K], T [K].
inline scalar FACChemistryModel::Kb(const SpeciesConst& sp, const scalar T) const
{
    const scalar R = chem_.Rcal;

    const scalar Tk = max(T, SMALL);

    const scalar lnKb =
        (-sp.Ab/Tk + sp.Bb*(log(Tk)-1.0) + sp.Db)/R;

    return exp(lnKb);
}

inline Foam::scalar Foam::FACChemistryModel::FmH
(
    const scalar mH,
    const scalar k0, const scalar k1, const scalar k2, const scalar k3,
    const scalar pH2, const scalar Kw, const scalar mA
) const
{
    (void)k2;
    // a = pH2^(1/3)
    // Use cbrt for speed and numerical stability
    const scalar a = cbrt(max(pH2, SMALL));

    const scalar mH2 = mH*mH;
    const scalar mH3 = mH2*mH;

    // F(mH) = LHS - RHS
    return 2.0*k0*a*mH3
         + (k1*a + 1.0)*mH2
         - mA*mH
         - (k3*a + Kw);
    // k2 is currently unused in this equation (kept in signature for consistency)
}

inline Foam::scalar Foam::FACChemistryModel::dFmH
(
    const scalar mH,
    const scalar k0, const scalar k1, const scalar k2, const scalar k3,
    const scalar pH2, const scalar Kw, const scalar mA
) const
{
    (void)k2; (void)k3; (void)Kw;

    const scalar a = cbrt(max(pH2, SMALL));

    // F'(mH) = derivative w.r.t mH
    return 6.0*k0*a*mH*mH
         + 2.0*(k1*a + 1.0)*mH
         - mA;

    // k2,k3,Kw do not contribute to derivative (except k3 in constant term)
}

void FACChemistryModel::compute
(
    const FACSelection& sel,
    const volScalarField& T,
    surfaceScalarField& CeqFace
) const
{
    // Initialize to zero (for non-selected faces)
    CeqFace = dimensionedScalar("zero", CeqFace.dimensions(), 0.0);

    const labelList& faces = sel.globalFaces();
    
    label nFail = 0;

    forAll(faces, i)
    {
        if (!sel.isBoundaryFace(i))
        {
            if (props_.skipInternalFaces) continue;

            FatalErrorInFunction
                << "Internal face in selection. FAC requires wall boundary faces."
                << exit(FatalError);
        }

        const label gf  = faces[i];
        const label own = mesh_.owner()[gf];

        const label pId = sel.patchId(i);
        const label lf  = sel.patchLocalFace(i);

        const scalar Tk   = T[own];
        const scalar pH2v = pH2Tab_->value(Tk);
        const scalar KwLog10 = KwTab_->value(Tk);          // table stores log10(Kw)
        const scalar Kwv     = pow(10.0, KwLog10);         // convert to linear Kw

        // species constants from chem_
        const scalar k0 = Kb(chem_.species[0], Tk);
        const scalar k1 = Kb(chem_.species[1], Tk);
        const scalar k2 = Kb(chem_.species[2], Tk);
        const scalar k3 = Kb(chem_.species[3], Tk);

        // Newton solve for mH
        /* scalar mH = 1e-8;
        mH = max(mH, newton_.xMin); */

        const scalar a = cbrt(max(pH2v, SMALL));
        const scalar rhs = max(k3*a + Kwv, SMALL);

        // dominant balance: 2*k0*a*mH^3 ~= rhs
        scalar mH = cbrt(rhs / max(2.0*k0*a, SMALL));

        // clamp to physical bounds
        mH = min(max(mH, newton_.xMin), newton_.xMax);

        auto F  = [&](const scalar x){ return FmH (x, k0,k1,k2,k3, pH2v,Kwv, mA_); };
        auto dF = [&](const scalar x){ return dFmH(x, k0,k1,k2,k3, pH2v,Kwv, mA_); };

        const bool ok = newtonSolve(mH, F, dF, newton_, false);
        if (!ok)
        {
            nFail++;
            if (nFail <= 5 && Pstream::master())
            {
                const scalar F0 = F(mH);
                const scalar d0 = dF(mH);

                WarningInFunction
                    << "Newton failed face i=" << i
                    << " own=" << own
                    << " T=" << Tk
                    << " pH2=" << pH2v
                    << " Kw=" << Kwv
                    << " k0=" << k0 << " k1=" << k1 << " k2=" << k2 << " k3=" << k3
                    << " mA=" << mA_
                    << " mH0=" << mH
                    << " F(mH0)=" << F0
                    << " dF(mH0)=" << d0
                    << nl << endl;
            }
            mH = max(mH, newton_.xMin);
        }

        //Find ferrous iron solubility Ceq in mol/kg
        const scalar Ceq = k0*sqr(mH)*a + (k1*mH)*a + k2*a + (k3/mH)*a;

        //k0.correctBoundaryConditions();
        CeqFace.boundaryFieldRef()[pId][lf] = Ceq;
    }
    reduce(nFail, sumOp<label>());

    if (nFail > 0 && Pstream::master())
    {
        WarningInFunction
            << "Newton mH solve failed on " << nFail
            << " selected faces. Using limited mH for those." << nl << endl;
    }
}

} // namespace Foam
