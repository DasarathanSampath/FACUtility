#include "FACChemistryTable.H"
#include "IFstream.H"
#include "dictionary.H"
#include "error.H"

namespace Foam
{

static dictionary readDict(const fileName& f)
{
    IFstream is(f);
    if (!is.good())
        FatalErrorInFunction << "Cannot open chemistry file: " << f << exit(FatalError);
    dictionary d(is);
    return d;
}

FACChemistryTable::FACChemistryTable(const fileName& chemistryDictFile)
:
    R_(1.987),
    Molor_mass_iron_(55.845),
    Density_iron_(7874),
    sp_(),
    Kw_(nullptr)
{
    dictionary d = readDict(chemistryDictFile);

    R_ = d.getOrDefault<scalar>("R_cal_per_molK", 1.987);
    Molor_mass_iron_ = d.getOrDefault<scalar>("Molor_mass_iron", 55.845);
    Density_iron_ = d.getOrDefault<scalar>("Density_iron", 7874);

    dictionary sc = d.subDict("speciesConstants");

    for (const word& k : wordList({"K0","K1","K2","K3"}))
    {
        if (!sc.found(k))
            FatalErrorInFunction << "Missing " << k << " in speciesConstants" << exit(FatalError);

        const List<scalar> vals(sc.lookup(k)); // reads "(A B D)"
        if (vals.size() != 3)
            FatalErrorInFunction << k << " must be (A B D)" << exit(FatalError);

        sp_.insert(k, SpeciesABC{vals[0], vals[1], vals[2]});
    }
}

void FACChemistryTable::setKwTable(const fileName& kwFile)
{
    Kw_.reset(new FACInterpolationTable(kwFile));
}

scalar FACChemistryTable::KfromABC(const SpeciesABC& abc, scalar T) const
{
    // lnKb = [-A/T + B*(ln(T)-1) + D]/R
    scalar lnK = (-abc.A/T + abc.B*(log(T) - 1.0) + abc.D)/R_;
    return exp(lnK);
}

void FACChemistryTable::constants
(
    scalar T,
    scalar& K0, scalar& K1, scalar& K2, scalar& K3,
    scalar& Kw
) const
{
    K0 = KfromABC(sp_["K0"], T);
    K1 = KfromABC(sp_["K1"], T);
    K2 = KfromABC(sp_["K2"], T);
    K3 = KfromABC(sp_["K3"], T);

    if (!Kw_.valid())
        FatalErrorInFunction << "Kw table not set" << exit(FatalError);

    Kw = Kw_->value(T);
}

} // namespace Foam
