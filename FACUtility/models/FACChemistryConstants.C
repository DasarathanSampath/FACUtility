#include "FACChemistryConstants.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "error.H"

namespace Foam
{

ChemistryConstants::ChemistryConstants()
:
    Rcal(1.987),
    MFe(55.845),
    rhoFe(7874),
    species()
{}

ChemistryConstants::ChemistryConstants(const fileName& f)
:
    ChemistryConstants()
{
    read(f);
}

static bool isCommentOrEmpty(const Foam::string& s)
{
    for (std::size_t i = 0; i < s.size(); ++i)
    {
        const char c = s[i];
        if (c == '#') return true;                 // comment line
        if (!std::isspace(static_cast<unsigned char>(c)))
            return false;                          // has content
    }
    return true;                                   // empty/whitespace
}

void ChemistryConstants::read(const fileName& f)
{
    IFstream is(f);
    if (!is.good())
    {
        FatalErrorInFunction << "Cannot open chemistry constants file: " << f
                             << exit(FatalError);
    }

    // Temporary storage
    DynamicList<SpeciesConst, 8> sp;

    bool inSpeciesBlock = false;

    while (!is.eof())
    {
        string line;
        is.getLine(line);

        if (line.empty() || isCommentOrEmpty(line)) continue;

        // Strip anything after '#'
        const label hashPos = line.find('#');
        if (hashPos != -1) line = line.substr(0, hashPos);

        // Detect block start/end
        if (line.find("speciesConstants") != string::npos) { inSpeciesBlock = true; continue; }
        if (inSpeciesBlock && line.find('{') != string::npos) continue;
        if (inSpeciesBlock && line.find('}') != string::npos) { inSpeciesBlock = false; continue; }

        IStringStream ls(line);

        if (!inSpeciesBlock)
        {
            word key;
            scalar val;
            if (ls.read(key) && ls.read(val))
            {
                if (key == "R_cal_per_molK")      Rcal = val;
                else if (key == "Molor_mass_iron") MFe = val;     // keep your spelling
                else if (key == "Density_iron")    rhoFe = val;
            }
        }
        else
        {
            scalar Ab, Bb, Db;
            if (ls.read(Ab) && ls.read(Bb) && ls.read(Db))
            {
                SpeciesConst sC;
                sC.Ab = Ab;
                sC.Bb = Bb;
                sC.Db = Db;
                sp.append(sC);
            }
        }
    }

    species.transfer(sp);

    if (species.size() != 4)
    {
        WarningInFunction
            << "Expected 4 species constants, found " << species.size()
            << ". Check file: " << f << nl << endl;
    }
}

} // namespace Foam
