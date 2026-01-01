#pragma once

#include "scalar.H"
#include "fileName.H"
#include "List.H"
#include "DynamicList.H"   // IMPORTANT for DynamicList<...>

namespace Foam
{

struct SpeciesConst
{
    scalar Ab{0};
    scalar Bb{0};
    scalar Db{0};
};

class ChemistryConstants
{
public:
    scalar Rcal{1.987};
    scalar MFe{55.845};
    scalar rhoFe{7.874};

    List<SpeciesConst> species;

    ChemistryConstants();
    explicit ChemistryConstants(const fileName& f);

    void read(const fileName& f);
};

} // namespace Foam
