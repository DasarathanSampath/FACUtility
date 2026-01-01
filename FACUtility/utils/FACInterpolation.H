#pragma once
#include "fileName.H"
#include "scalarList.H"

namespace Foam
{
class FACInterpolationTable
{
    scalarList x_, y_;

public:
    explicit FACInterpolationTable(const fileName& f);
    scalar value(scalar xq) const; // linear clamp+interp
};
}
