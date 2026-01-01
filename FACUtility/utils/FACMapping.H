#pragma once

#include "fvMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "FACSelection.H"
#include "error.H"

namespace Foam
{

/**
 * Map selected boundary faces (from FACSelection) to owner cells.
 *
 * Modes:
 *  - max:   cell = max over selected faces connected to that owner cell
 *  - mean:  cell = average over selected faces connected to that owner cell
 *
 * Notes:
 *  - Only maps boundary faces. Internal faces are skipped (or error if desired).
 *  - Only maps for selected faces. Other cells remain unchanged unless you zero first.
 */
enum class MapMode { MAX, MEAN };

inline void mapSurfaceToVolOwner
(
    const fvMesh& mesh,
    const FACSelection& sel,
    const surfaceScalarField& sField,
    volScalarField& vField,
    const MapMode mode,
    const bool zeroBefore = true,
    const bool skipInternal = true
)
{
    if (zeroBefore)
    {
        vField = dimensionedScalar("zero", vField.dimensions(), 0.0);
    }

    // For mean mode, track counts
    labelList count;
    if (mode == MapMode::MEAN)
    {
        count.setSize(mesh.nCells(), 0);
    }

    const labelList& faces = sel.globalFaces();

    forAll(faces, i)
    {
        const label gf = faces[i];

        if (!sel.isBoundaryFace(i))
        {
            if (skipInternal) continue;
            FatalErrorInFunction
                << "Internal face encountered in selection during mapping. "
                << "Set skipInternalFaces=true or use a wall-only faceZone."
                << exit(FatalError);
        }

        const label pId = sel.patchId(i);
        const label lf  = sel.patchLocalFace(i);

        const scalar val = sField.boundaryField()[pId][lf];
        const label owner = mesh.faceOwner()[gf];

        if (mode == MapMode::MAX)
        {
            vField[owner] = max(vField[owner], val);
        }
        else // MEAN
        {
            vField[owner] += val;
            count[owner] += 1;
        }
    }

    if (mode == MapMode::MEAN)
    {
        forAll(count, c)
        {
            if (count[c] > 0)
            {
                vField[c] /= scalar(count[c]);
            }
        }
    }

    vField.correctBoundaryConditions(); // usually harmless for volScalarField post fields
}

/**
 * Map selected boundary faces (from FACSelection) to a *volScalarField patch*.
 *
 * Purpose:
 *  - ParaView often won’t show surfaceScalarField nicely.
 *  - This creates a normal volScalarField (e.g., TLpatch) whose patch values
 *    match the face field on the selected faces.
 *
 * Behavior:
 *  - internalField is set to internalFill (default 0)
 *  - all boundary patches are set to internalFill first
 *  - then ONLY selected faces are filled from sField
 *
 * Notes:
 *  - Do NOT call correctBoundaryConditions() here by default, because a patch type
 *    could overwrite what we just assigned. We want exact values for visualization.
 */
inline void mapSurfaceToVolPatch
(
    const fvMesh& mesh,
    const FACSelection& sel,
    const surfaceScalarField& sField,
    volScalarField& vField,
    const bool zeroBefore = true,
    const bool skipInternal = true
)
{
    if (zeroBefore)
    {
        vField = dimensionedScalar("zero", vField.dimensions(), 0.0);
    }

    // 3) overwrite only selected boundary faces
    const labelList& faces = sel.globalFaces();

    forAll(faces, i)
    {
        if (!sel.isBoundaryFace(i))
        {
            if (skipInternal) continue;
            FatalErrorInFunction
                << "Internal face encountered in selection during patch mapping."
                << exit(FatalError);
        }

        const label pId = sel.patchId(i);
        const label lf  = sel.patchLocalFace(i);

        vField.boundaryFieldRef()[pId][lf] = sField.boundaryField()[pId][lf];
    }

    // Intentionally no vField.correctBoundaryConditions();
}

} // namespace Foam
