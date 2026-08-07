#pragma once

#include "fvMesh.H"
#include "surfaceFields.H"
#include "FACSelection.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"

namespace Foam
{

inline fileName makeDatPath
(
    const Time& runTime,
    const word& tag              // e.g. patch name or faceZone name
)
{
    // postProcessing/FACUtility/<time>/FAC_<tag>[_procN].dat
    fileName dir = runTime.path()/"postProcessing"/"FACUtility"/runTime.timeName();
    mkDir(dir);

    return dir / ("FAC_" + tag + ".dat");
}

/**
 * Write per-selected-face values with spatial locations.
 * Writes boundary faces only (internal faces skipped).
 *
 * Columns:
 * x y z UTau km Ceq TL temperature
 */
inline void writeFACDat
(
    const fvMesh& mesh,
    const FACSelection& sel,
    const Time& runTime,
    const word& tag,
    const surfaceScalarField& UTauFace,
    const surfaceScalarField& kmFace,
    const surfaceScalarField& CeqFace,
    const surfaceScalarField& TLFace,
    const volScalarField& T,
    const volScalarField& yPlus,
    const bool skipInternal = true
)
{
    // Avoid multiple writers in parallel
    if (Pstream::parRun() && !Pstream::master())
    {
        return;
    }
    // Choose path + open file (overwrite each time-step)
    const fileName outFile = makeDatPath(runTime, tag);

    OFstream os
    (
        outFile,
        IOstream::ASCII,
        IOstream::UNCOMPRESSED,
        IOstreamOption::NON_APPEND
    );

    // Header
    os  << "# FACUtility face data\n";
    os  << "# time = " << runTime.timeName() << "\n";
    os  << "# columns: x y z UTau km Ceq TL T yPlus\n";

    const labelList& faces = sel.globalFaces();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(faces, i)
    {
        if (!sel.isBoundaryFace(i))
        {
            if (skipInternal) continue;
            // internal faces are not meaningful for wall FAC
            continue;
        }

        const label pId = sel.patchId(i);
        const label lf  = sel.patchLocalFace(i);

        // Face center on this patch face
        // (polyPatch::faceCentres is boundary-face centres in global coords)
        const polyPatch& pp = pbm[pId];
        const point& c = pp.faceCentres()[lf];

        const scalar ut  = UTauFace.boundaryField()[pId][lf];
        const scalar km  = kmFace.boundaryField()[pId][lf];
        const scalar ceq = CeqFace.boundaryField()[pId][lf];
        const scalar tl  = TLFace.boundaryField()[pId][lf];
        const scalar temperature  = T.boundaryField()[pId][lf];
        const scalar yPlusValue = yPlus.boundaryField()[pId][lf];

        os  << c.x() << ' ' << c.y() << ' ' << c.z() << ' '
            << ut    << ' ' << km    << ' ' << ceq   << ' ' << tl << ' ' 
            << temperature << ' ' << yPlusValue
            << '\n';
    }
}

} // namespace Foam
