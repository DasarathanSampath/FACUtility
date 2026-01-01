#include "FACSelection.H"
#include "FACProperties.H"
#include "faceZoneMesh.H"
#include "polyBoundaryMesh.H"
#include "error.H"

namespace Foam
{

FACSelection::FACSelection(const fvMesh& mesh, const FACProperties& props)
:
    mesh_(mesh),
    props_(props),
    globalFaces_(),
    isBoundary_(),
    patchId_(),
    patchLocalFace_()
{
    if (props_.selectionMode == "patch")
        buildFromPatch();
    else if (props_.selectionMode == "faceZone")
        buildFromFaceZone();
    else
        FatalErrorInFunction << "selectionMode must be patch or faceZone" << exit(FatalError);

    buildBoundaryMap();
}

void FACSelection::buildFromPatch()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    label pId = pbm.findPatchID(props_.patchName);

    if (pId < 0)
        FatalErrorInFunction << "Patch not found: " << props_.patchName << exit(FatalError);

    const polyPatch& pp = pbm[pId];
    globalFaces_.setSize(pp.size());

    label start = pp.start();
    forAll(globalFaces_, i)
        globalFaces_[i] = start + i;
}

void FACSelection::buildFromFaceZone()
{
    const faceZoneMesh& fzm = mesh_.faceZones();
    label zId = fzm.findZoneID(props_.faceZoneName);

    if (zId < 0)
        FatalErrorInFunction << "faceZone not found: " << props_.faceZoneName << exit(FatalError);

    globalFaces_ = fzm[zId]; // global face labels
}
/*-------*\
For each selected face (from patch or faceZone),
is it a boundary face?
This is important because faceZones may include internal faces.
If yes:
which patch does it belong to?
what is its local face index inside that patch?
\*-------*/

void FACSelection::buildBoundaryMap()
{
    isBoundary_.setSize(globalFaces_.size(), false); // Is this selected face a boundary face?
    patchId_.setSize(globalFaces_.size(), -1); // Which boundary patch it belongs to
    patchLocalFace_.setSize(globalFaces_.size(), -1); // Index inside that patch

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh(); //Access the boundary mesh

    forAll(globalFaces_, i)
    {
        label gf = globalFaces_[i];

        if (mesh_.isInternalFace(gf)) //Skip internal faces
            continue;

        label pId = pbm.whichPatch(gf); //Find which patch owns this face
        if (pId < 0) continue;

        // Convert global face index → local patch face index
        const polyPatch& pp = pbm[pId];
        isBoundary_[i] = true;
        patchId_[i] = pId;
        patchLocalFace_[i] = gf - pp.start();
    }
}

} // namespace Foam
