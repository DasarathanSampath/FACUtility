// Required fields in time directory
facProps.UName = FACDict.getOrDefault<word>("UName", "U");
facProps.TName = FACDict.getOrDefault<word>("TName", "T");
facProps.pName = FACDict.getOrDefault<word>("pName", "p");

// ---- Vicosity details ------
facProps.rhoName  = FACDict.getOrDefault<word>("rhoName", "rho");
facProps.muEffName = FACDict.getOrDefault<word>("muEffName", "muEff");
facProps.useThermoForRho = FACDict.getOrDefault<Switch>("useThermoForRho", true);

// ---- Hydraulic diameter / bulk flow velocity ----
// facProps.charLength = readScalar(FACDict.lookup("charLength"));
facProps.charLength = FACDict.getOrDefault<scalar>("charLength", 0.005);
facProps.UmeanConstant = FACDict.getOrDefault<scalar>("UmeanConstant", 0.005);

// ---- Ferrous diffusivity D: constant or table(T) ---- 
facProps.diffusivityMode = FACDict.getOrDefault<word>("diffusivityMode", "constant");
facProps.Dconstant = FACDict.getOrDefault<scalar>("Dconstant", 1e-8);
facProps.DtableFile = FACDict.getOrDefault<fileName>("DtableFile", "");

// ---- Chemistry inputs ----
facProps.mA = FACDict.getOrDefault<scalar>("mA", 7.94e-3);
facProps.chromiumFactor = FACDict.getOrDefault<scalar>("chromiumFactor", 0.947);
facProps.corrosionHours = FACDict.getOrDefault<scalar>("corrosionHours", 10000);
facProps.chemistryTableFile = FACDict.getOrDefault<fileName>("chemistryTableFile", "");
facProps.KwTableFile = FACDict.getOrDefault<fileName>("KwTableFile", "");
facProps.pH2TableFile = FACDict.getOrDefault<fileName>("pH2TableFile", "");

// ---- Output controls ----
facProps.writeCSV = FACDict.getOrDefault<Switch>("writeCSV", true);
facProps.csvFile = FACDict.getOrDefault<fileName>("csvFile", "postProcessing/FACUtility/FACsummary.csv");

// ---- Other defaults ----
facProps.skipInternalFaces = FACDict.getOrDefault<Switch>("skipInternalFaces", true);
facProps.writeVolMapping = FACDict.getOrDefault<Switch>("writeVolMapping", true);

// ---- Region selection: patch or cellZone ----
facProps.selectionMode = FACDict.getOrDefault<word>("selectionMode", "patch");
facProps.patchName = FACDict.getOrDefault<word>("patchName", "wall");
facProps.faceZoneName = FACDict.getOrDefault<word>("faceZoneName", "");

if (FACDict.found("newton"))
{
    const dictionary& nd = FACDict.subDict("newton");

    nd.readIfPresent("maxIter",     facProps.newton.maxIter);
    nd.readIfPresent("absTol",      facProps.newton.absTol);
    nd.readIfPresent("relTol",      facProps.newton.relTol);
    nd.readIfPresent("minStep",     facProps.newton.minStep);
    nd.readIfPresent("maxDampIter", facProps.newton.maxDampIter);
    nd.readIfPresent("xMin",        facProps.newton.xMin);
    nd.readIfPresent("xMax",        facProps.newton.xMax);
}
