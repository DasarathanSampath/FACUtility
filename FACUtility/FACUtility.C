// OpenFOAM Native headers
#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "fluidThermo.H"
#include "basicThermo.H"

// FAC input data processing headers
#include "FACProperties.H"
#include "FACSelection.H"

// FAC prediction model headers
#include "FACErosionModel.H"
#include "FACChemistryModel.H"
#include "FACThinningModel.H"

// FAC output data processing headers
#include "FACMapping.H"
#include "FACCSV.H"
#include "FACDatWriter.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"

    // read FAC input data from systems/FACProperties
    IOdictionary FACDict
    (
        IOobject("FACProperties", runTime.system(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    FACProperties facProps;
    #include "readFACProperties.H"

    // create folder for postProcessing outputs
    mkDir(runTime.path()/"postProcessing"/"FACUtility");


    forAll(timeDirs, ti)
    {
        runTime.setTime(timeDirs[ti], ti);
        Info<< "\nTime = " << runTime.timeName() << nl << endl;
        mesh.readUpdate();

        // create FAC fields
        #include "createFields.H"

        // Selection of patch/faceZone
        FACSelection sel(mesh, facProps);
        Info << "\n Cell Selection done\n" << endl;

        // FAC Models
        FACErosionModel erosion(mesh, facProps);
        FACChemistryModel chemistry(mesh, facProps);
        FACThinningModel thinning(mesh, facProps);
        Info << "\n The models execution done\n" << endl;

        // Compute FAC surface fields
        erosion.compute(sel, U, T, nuEff, UTauFace, kmFace);
        chemistry.compute(sel, T, CeqFace);
        thinning.compute(sel, rho, kmFace, CeqFace, TLFace);
        Info << "\n The models Compute done\n" << endl;

        // write FAC field data
        #include "writeOutputs.H"
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}
