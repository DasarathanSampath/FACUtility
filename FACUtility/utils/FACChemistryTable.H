#pragma once
#include "fileName.H"
#include "HashTable.H"
#include "FACInterpolation.H"

namespace Foam
{

struct SpeciesABC { scalar A, B, D; };

class FACChemistryTable
{
    scalar R_; // cal/mol-K
    scalar Molor_mass_iron_; // g/mol
    scalar Density_iron_; // g/cm³
    HashTable<SpeciesABC> sp_; // K0..K3
    autoPtr<FACInterpolationTable> Kw_; // required

public:
    explicit FACChemistryTable(const fileName& chemistryDictFile);

    void setKwTable(const fileName& kwFile);

    void constants
    (
        scalar T,
        scalar& K0, scalar& K1, scalar& K2, scalar& K3,
        scalar& Kw
    ) const;

private:
    scalar KfromABC(const SpeciesABC& abc, scalar T) const;
};

}
