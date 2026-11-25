#!/usr/bin/env python3
"""
Make small owner/boundary files for one wall patch, for the exact
face-by-face check in verify_cases.py.

The full constant/polyMesh/owner file lists the owner cell of every face in
the mesh, so it can be very large. verify_cases.py only needs the owner cells
of the patch that FACUtility writes (patchName in system/FACProperties).
This script copies only those entries. The cell numbers are not changed, so
the small files work with the full T file of the same case.

Usage (run on the reconstructed case, not on processor* folders):

    python3 extract_patch_owner.py <caseDir> <patchName> <outDir>

    example:
    python3 extract_patch_owner.py ~/run/case4 walls_elbow nFail_verification/case4

It writes <outDir>/owner and <outDir>/boundary. Then run verify_cases.py as
before. It now reports "owner-cell T (mesh owner/boundary files)" and the
maximum face-by-face relative difference for each case.

Both ascii and binary owner files are read (32- or 64-bit labels).
"""

import os
import re
import struct
import sys


def strip_comments(t):
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    return re.sub(r"//[^\n]*", "", t)


def patch_range(boundary_file, patch):
    t = strip_comments(open(boundary_file, errors="replace").read())
    m = re.search(r"(?<![\w.])" + re.escape(patch) + r"\s*\{([^}]*)\}", t)
    if not m:
        names = re.findall(r"(\w+)\s*\{[^}]*nFaces", t)
        sys.exit(f"patch '{patch}' not found in {boundary_file}. Patches: {', '.join(names)}")
    start = int(re.search(r"startFace\s+(\d+)", m.group(1)).group(1))
    n = int(re.search(r"nFaces\s+(\d+)", m.group(1)).group(1))
    return start, n


def read_owner_slice(owner_file, start, n):
    with open(owner_file, "rb") as f:
        head = f.read(4096)
        hdr_end = head.find(b"}")
        if hdr_end < 0:
            sys.exit(f"{owner_file}: FoamFile header not found")
        header = head[:hdr_end].decode(errors="replace")
        binary = re.search(r"format\s+binary", header) is not None
        lsize = 8 if re.search(r"label\s*=\s*64", header) else 4

        # after the header: comment lines, then "<count>" and "(" (OpenFOAM
        # writes them on separate lines, or as "<count>(" for short lists)
        f.seek(hdr_end + 1)
        buf = f.read(65536)
        m = re.search(rb"(?m)^[ \t]*(\d+)[ \t]*\r?\n?[ \t]*\(", buf)
        if not m:
            sys.exit(f"{owner_file}: list start not found")
        total = int(m.group(1))
        if start + n > total:
            sys.exit(f"{owner_file}: has {total} entries, patch needs {start + n}")
        list_start = hdr_end + 1 + m.end()

        if binary:
            f.seek(list_start + start * lsize)
            data = f.read(n * lsize)
            fmt = "<%d%s" % (n, "q" if lsize == 8 else "i")
            return list(struct.unpack(fmt, data)), total

        # ascii: one label per entry, separated by white space
        f.seek(list_start)
        out, idx, rest = [], 0, b""
        while len(out) < n:
            chunk = f.read(1 << 22)
            if not chunk:
                break
            parts = (rest + chunk).split()
            rest = parts.pop() if chunk[-1:] not in b" \t\r\n" else b""
            for p in parts:
                if p == b")":
                    break
                if start <= idx < start + n:
                    out.append(int(p))
                idx += 1
                if len(out) == n:
                    break
        if len(out) < n and rest and rest != b")":
            out.append(int(rest))
        if len(out) != n:
            sys.exit(f"{owner_file}: read {len(out)} of {n} patch entries")
        return out, total


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    case, patch, outdir = sys.argv[1:]
    poly = os.path.join(case, "constant", "polyMesh")
    if not os.path.isdir(poly):
        sys.exit(f"{poly} not found (use the reconstructed case)")
    start, n = patch_range(os.path.join(poly, "boundary"), patch)
    own, total = read_owner_slice(os.path.join(poly, "owner"), start, n)

    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "owner"), "w") as f:
        f.write("FoamFile\n{\n    format      ascii;\n    class       labelList;\n"
                "    object      owner;\n}\n")
        f.write(f"// owner cells of patch '{patch}' only (faces {start}-{start + n - 1} "
                f"of {total}), made by extract_patch_owner.py\n\n")
        f.write(f"{n}\n(\n" + "\n".join(map(str, own)) + "\n)\n")
    with open(os.path.join(outdir, "boundary"), "w") as f:
        f.write("FoamFile\n{\n    format      ascii;\n    class       polyBoundaryMesh;\n"
                "    object      boundary;\n}\n\n1\n(\n")
        f.write(f"    {patch}\n    {{\n        type            wall;\n"
                f"        nFaces          {n};\n        startFace       0;\n    }}\n)\n")
    print(f"{patch}: {n} faces (startFace {start} in a mesh of {total} faces), "
          f"owner cells {min(own)}-{max(own)} -> {outdir}/owner, {outdir}/boundary")


if __name__ == "__main__":
    main()
