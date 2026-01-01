#include "FACInterpolation.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "error.H"
#include "DynamicList.H"

namespace Foam
{

FACInterpolationTable::FACInterpolationTable(const fileName& f)
:
    x_(),
    y_()
{
    IFstream is(f);
    if (!is.good())
    {
        FatalErrorInFunction << "Cannot open table: " << f << exit(FatalError);
    }

    DynamicList<scalar, 16> X, Y;   // FIX for OpenFOAM-2212

    while (!is.eof())
    {
        string line;
        is.getLine(line);
        if (line.empty() || line[0] == '#') continue;

        IStringStream ls(line);
        scalar a, b;
        if (ls.read(a) && ls.read(b))
        {
            X.append(a);
            Y.append(b);
        }
    }

    if (X.size() < 2)
    {
        FatalErrorInFunction << "Table needs >=2 points: " << f << exit(FatalError);
    }

    x_.transfer(X);
    y_.transfer(Y);
}

scalar FACInterpolationTable::value(scalar xq) const
{
    if (xq <= x_.first()) return y_.first();
    if (xq >= x_.last())  return y_.last();

    for (label i=0; i<x_.size()-1; ++i)
    {
        if (xq >= x_[i] && xq <= x_[i+1])
        {
            scalar t = (xq - x_[i])/(x_[i+1] - x_[i]);
            return (1.0 - t)*y_[i] + t*y_[i+1];
        }
    }
    return y_.last();
}

} // namespace Foam
