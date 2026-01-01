#pragma once
#include "fvMesh.H"
#include "labelList.H"
#include "boolList.H"

namespace Foam
{
struct FACProperties;

class FACSelection
{
    const fvMesh& mesh_;
    const FACProperties& props_;

    labelList globalFaces_;     // global face labels
    boolList isBoundary_;       // per selected face
    labelList patchId_;         // patch id if boundary else -1
    labelList patchLocalFace_;  // local face index in patch if boundary else -1

    void buildFromPatch();
    void buildFromFaceZone();
    void buildBoundaryMap();

public:
    FACSelection(const fvMesh& mesh, const FACProperties& props);

    const labelList& globalFaces() const { return globalFaces_; }
    bool isBoundaryFace(label i) const { return isBoundary_[i]; }
    label patchId(label i) const { return patchId_[i]; }
    label patchLocalFace(label i) const { return patchLocalFace_[i]; }

    label size() const { return globalFaces_.size(); }
};

} // namespace Foam
