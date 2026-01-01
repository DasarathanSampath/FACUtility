#pragma once

#include "fvMesh.H"
#include "FACSelection.H"
#include "surfaceFields.H"
#include "OFstream.H"
#include "fileName.H"
#include "OSspecific.H"   // for isFile / mkDir
#include "PstreamReduceOps.H" // if you later add parallel reductions (optional)

namespace Foam
{

inline void appendFACCSV
(
    const fvMesh& mesh,
    const FACSelection& sel,
    const fileName& csvFile,
    const word& timeName,
    const surfaceScalarField& UTauFace,
    const surfaceScalarField& kmFace,
    const surfaceScalarField& CeqFace,
    const surfaceScalarField& TLFace,
    const bool skipInternal = true
)
{
    // Ensure directory exists
    mkDir(csvFile.path());

    const bool newFile = !isFile(csvFile);

    OFstream os
    (
        csvFile,
        IOstream::ASCII,
        IOstream::UNCOMPRESSED,
        IOstream::APPEND
    );


    if (newFile)
    {
        os  << "time,"
            << "UTau_mean,UTau_max,"
            << "km_mean,km_max,"
            << "Ceq_mean,Ceq_max,"
            << "TL_mean,TL_max\n";
    }

    scalar sumUT=0, sumKm=0, sumCeq=0, sumTL=0;
    scalar maxUT=0, maxKm=0, maxCeq=0, maxTL=0;
    label n = 0;

    const labelList& faces = sel.globalFaces();

    forAll(faces, i)
    {
        if (!sel.isBoundaryFace(i))
        {
            if (skipInternal) continue;
            // If not skipping, you can treat internal as invalid for FAC
            continue;
        }

        const label pId = sel.patchId(i);
        const label lf  = sel.patchLocalFace(i);

        const scalar ut  = UTauFace.boundaryField()[pId][lf];
        const scalar km  = kmFace.boundaryField()[pId][lf];
        const scalar ceq = CeqFace.boundaryField()[pId][lf];
        const scalar tl  = TLFace.boundaryField()[pId][lf];

        sumUT += ut;  sumKm += km;  sumCeq += ceq;  sumTL += tl;
        maxUT = max(maxUT, ut);
        maxKm = max(maxKm, km);
        maxCeq = max(maxCeq, ceq);
        maxTL = max(maxTL, tl);
        n++;
    }

    const scalar UTmean  = (n > 0 ? sumUT / scalar(n) : 0.0);
    const scalar Kmmean  = (n > 0 ? sumKm / scalar(n) : 0.0);
    const scalar Ceqmean = (n > 0 ? sumCeq / scalar(n) : 0.0);
    const scalar TLmean  = (n > 0 ? sumTL / scalar(n) : 0.0);

    os  << timeName << ","
        << UTmean  << "," << maxUT  << ","
        << Kmmean  << "," << maxKm  << ","
        << Ceqmean << "," << maxCeq << ","
        << TLmean  << "," << maxTL  << "\n";
}

} // namespace Foam
