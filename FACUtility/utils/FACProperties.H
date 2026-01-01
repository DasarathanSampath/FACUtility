#pragma once
#include "word.H"
#include "fileName.H"
#include "Switch.H"

#include "FACNewtonRaphson.H"

namespace Foam
{

struct FACProperties
{
    // Required fields in time directory
    word UName {"U"};
    word TName {"T"};
    word pName {"p"};
    word rhoName  {"rho"};
    word muEffName {"muEff"};
    Switch useThermoForRho {true};

    // ---- Hydraulic diameter / bulk flow velocity ----
    scalar charLength {0.005};        // [m] reference diameter / hydraulic length
    scalar UmeanConstant {4.5};        // [m/s] bulk velocity

    // ---- Ferrous diffusivity D: constant or table(T) ---- 
    word diffusivityMode {"constant"}; // constant | table
    scalar Dconstant {1e-8};           // [m2/s]
    fileName DtableFile {""};

    // ---- Chemistry inputs ----
    fileName chemistryTableFile {""};
    fileName KwTableFile {""};
    fileName pH2TableFile {""};        // optional
    scalar mA {7.94e-3};               // molality input
    scalar chromiumFactor {0.947};     // (1-beta)
    scalar corrosionHours {10000};

    // Newton-Raphson
    NewtonControls newton;

    // ---- Output controls ----
    Switch writeCSV {true};
    fileName csvFile {"postProcessing/FACUtility/FACsummary.csv"};

    // ---- Other defaults ----
    Switch writeVolMapping {true};     // also create volScalar mapped fields
    Switch skipInternalFaces {true};   // for faceZone: skip internal faces

    // Selection
    word selectionMode {"patch"};     // patch | faceZone
    word patchName {"wall"};
    word faceZoneName {""};
};

} // namespace Foam
