UTauFace.write();
kmFace.write();
CeqFace.write();
TLFace.write();

if (facProps.writeVolMapping)
{
    // Map selected wall faces -> owner cells
    mapSurfaceToVolOwner(mesh, sel, UTauFace, UTau, Foam::MapMode::MAX, true, facProps.skipInternalFaces);
    mapSurfaceToVolOwner(mesh, sel, kmFace, km, Foam::MapMode::MAX, true, facProps.skipInternalFaces);
    mapSurfaceToVolOwner(mesh, sel, CeqFace, Ceq, Foam::MapMode::MAX, true, facProps.skipInternalFaces);
    mapSurfaceToVolOwner(mesh, sel, TLFace, TL,   Foam::MapMode::MAX, true, facProps.skipInternalFaces);

    UTau.write();
    km.write();
    Ceq.write();
    TL.write();

    // Following Patch writes are only for visualization
    mapSurfaceToVolPatch(mesh, sel, UTauFace, UTauPatch, true, facProps.skipInternalFaces);
    mapSurfaceToVolPatch(mesh, sel, kmFace, kmPatch, true, facProps.skipInternalFaces);
    mapSurfaceToVolPatch(mesh, sel, CeqFace, CeqPatch, true, facProps.skipInternalFaces);
    mapSurfaceToVolPatch(mesh, sel, TLFace, TLPatch, true, facProps.skipInternalFaces);

    UTauPatch.write();
    kmPatch.write();
    CeqPatch.write();
    TLPatch.write();
}

if (facProps.writeCSV)
{
    appendFACCSV
    (
        mesh,
        sel,
        facProps.csvFile,
        runTime.timeName(),
        UTauFace,
        kmFace,
        CeqFace,
        TLFace,
        facProps.skipInternalFaces
    );
}

// Write .dat with full per-face values + xyz
word tag = (facProps.selectionMode == "patch" ? facProps.patchName : facProps.faceZoneName);

writeFACDat
(
    mesh,
    sel,
    runTime,
    tag,
    UTauFace,
    kmFace,
    CeqFace,
    TLFace,
    T,
    yPlus,
    facProps.skipInternalFaces
);
