#!/usr/bin/env python3
"""
Face-by-face check of CeqFace from the validation runs (Cases 1-4) against the
independent solution of Eqs. (6)-(9) in verify_chemistry.py, plus a check for
faces carrying the Newton-Raphson fallback value (nFail).

Folder layout (one subfolder per case; only small files are needed):

    nFail_verification/
      case1/CeqFace          required  (from the time directory FACUtility wrote)
      case1/T                required  (same time directory)
      case1/owner            optional  (owner cells of the patch faces, made by
                                        extract_patch_owner.py, or the full
                                        constant/polyMesh/owner in ascii)
      case1/boundary         optional  (matching boundary file)
      case1/FACProperties    required  (system/FACProperties: mA, patchName)
      Kw_vs_T.dat, pH2_vs_T.dat   optional, in the case folder or the root folder
                             (if absent, Listings 3 and 4 are used)
      case2/ ...

Run from that folder:
    python3 verify_cases.py            (or: python3 verify_cases.py --root <folder>)

FACUtility evaluates each face at its owner-cell temperature (Listing 8,
line 22). With owner and boundary present the check is exact, face by face.
Without them, every face is checked against the range of Ceq spanned by the
minimum and maximum cell temperature in the domain (a bounding check, which is
tight when T is nearly uniform, as in Cases 1-3).

Fallback detection: a failed solve clamps m(H+) to newton.xMin (1e-20), which
makes Ceq of order 1e4 mol/kg. Any face with Ceq > 1e-2 mol/kg, or non-finite
or non-positive, is counted as a fallback face.
"""

import argparse
import math
import os
import re
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify_chemistry as vc  # noqa: E402
from collect_ceq import patch_values, strip_comments  # noqa: E402

mp.mp.dps = 30


def read_props(path):
    t = strip_comments(open(path).read())
    def get(key):
        m = re.search(r"^\s*" + key + r"\s+([^;]+);", t, flags=re.M)
        return m.group(1).strip().strip('"') if m else None
    mA = get("mA")
    if mA is None:
        raise ValueError(f"{path}: mA not found")
    m = re.search(r"newton\s*\{([^}]*)\}", t)
    tol = re.search(r"absTol\s+([^;]+);", m.group(1)).group(1).strip() if m and "absTol" in m.group(1) else "?"
    return float(mA.replace(" ", "")), get("patchName"), get("selectionMode"), tol


def read_table(path):
    rows = []
    for line in open(path):
        line = line.split("#", 1)[0].split("//", 1)[0].strip()
        nums = line.replace(",", " ").split()
        if len(nums) >= 2:
            try:
                rows.append((float(nums[0]), float(nums[1])))
            except ValueError:
                pass
    if len(rows) < 2:
        raise ValueError(f"{path}: fewer than two numeric rows")
    return rows


def internal_values(path):
    """internalField of a volScalarField as a numpy array (ascii format)."""
    import numpy as np
    with open(path) as f:
        text = f.read()
    i = text.index("internalField")
    head = text[i:i + 200]
    m = re.match(r"internalField\s+uniform\s+([^;\s]+)", head)
    if m:
        return np.array([float(m.group(1))])
    m = re.match(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", text[i:i + 300])
    if not m:
        raise ValueError(f"{path}: cannot parse internalField (binary format? run foamFormatConvert)")
    n = int(m.group(1))
    start = i + m.end()
    end = text.index(")", start)
    vals = np.array(text[start:end].split(), dtype=float)
    if len(vals) != n:
        raise ValueError(f"{path}: expected {n} values, read {len(vals)}")
    return vals


def owner_cells(d, patch, nfaces):
    """Owner cells of the patch faces from <case>/owner and <case>/boundary, if present."""
    ofile, bfile = os.path.join(d, "owner"), os.path.join(d, "boundary")
    if not (os.path.isfile(ofile) and os.path.isfile(bfile)):
        return None
    b = strip_comments(open(bfile).read())
    m = re.search(r"\b" + re.escape(patch) + r"\s*\{([^}]*)\}", b)
    if not m:
        raise ValueError(f"boundary: patch '{patch}' not found")
    start = int(re.search(r"startFace\s+(\d+)", m.group(1)).group(1))
    n = int(re.search(r"nFaces\s+(\d+)", m.group(1)).group(1))
    if n != nfaces:
        raise ValueError(f"boundary says {n} faces on '{patch}', CeqFace has {nfaces}")
    with open(ofile, "rb") as f:
        raw = f.read()
    if b"format      binary" in raw[:2000] or b"format binary" in raw[:2000]:
        raise ValueError("owner is in binary format; run foamFormatConvert or copy an ascii one")
    text = strip_comments(raw.decode())
    text = text[text.index("}", text.index("FoamFile")) + 1:] if "FoamFile" in text else text
    m = re.search(r"(\d+)\s*\(", text)
    body = text[m.end():text.index(")", m.end())].split()
    return [int(x) for x in body[start:start + n]]


def sig_digits(tok):
    mant = re.split(r"[eE]", tok)[0].lstrip("+-").replace(".", "").lstrip("0")
    return len(mant)


def check_case(d, patch_override, tol_user):
    name = os.path.basename(d.rstrip("/"))
    mA, patch, mode, absTol = read_props(os.path.join(d, "FACProperties"))
    patch = patch_override or patch
    if mode and mode != "patch" and not patch_override:
        raise ValueError(f"selectionMode is '{mode}'; pass --patch with the patch to read")
    for fname, table in (("Kw_vs_T.dat", vc.KW_TABLE), ("pH2_vs_T.dat", vc.PH2_TABLE)):
        for cand in (os.path.join(d, fname), os.path.join(os.path.dirname(os.path.abspath(d)), fname)):
            if os.path.isfile(cand):
                table[:] = read_table(cand)
                break

    ceq, toks = patch_values(os.path.join(d, "CeqFace"), patch)
    digits = max(sig_digits(t) for t in toks)
    # rounding of the written Ceq (half a unit in the last digit) plus the rounding
    # of the written T (0.001 K gives up to ~4e-6 in Ceq): about twice the first term
    tol = tol_user or (1e-6 if digits >= 9 else 2 * 0.5 * 10 ** (1 - digits) * 1.01)

    bad = [i for i, c in enumerate(ceq) if not math.isfinite(c) or c <= 0 or c > 1e-2]

    # FACUtility evaluates each face at its owner-cell temperature, T[own]
    # (Listing 8, line 22), so wall-face temperatures come from the internal field.
    Tsrc, Tl, Trange = None, None, None
    tp = os.path.join(d, "T_patch.txt")
    if os.path.isfile(os.path.join(d, "T")):
        Tint = internal_values(os.path.join(d, "T"))
        own = owner_cells(d, patch, len(ceq))
        if own is not None:
            Tl = [float(Tint[c]) for c in own] if len(Tint) > 1 else [float(Tint[0])] * len(ceq)
            Tsrc = "owner-cell T (mesh owner/boundary files)"
        else:
            Trange = (float(Tint.min()), float(Tint.max()))
            Tsrc = "range of T over all cells (no owner/boundary files: bounding check)"
    elif os.path.isfile(tp):
        Tl = [float(open(tp).read().split()[0])] * len(ceq)
        Tsrc = "T_patch.txt (single value)"
    else:
        raise ValueError("no T file (or T_patch.txt) in the case folder")

    print(f"\n== {name}: patch '{patch}', {len(ceq)} faces, mA = {mA:g} mol/kg, "
          f"{digits} significant digits written -> tolerance {tol:.1e}")
    if Tl is not None:
        print(f"   wall T [{Tsrc}]: {min(Tl):.3f} - {max(Tl):.3f} K")
    else:
        print(f"   T [{Tsrc}]: {Trange[0]:.3f} - {Trange[1]:.3f} K")
    print(f"   newton.absTol in this FACProperties: {absTol}")
    print(f"   CeqFace: {min(ceq):.6e} - {max(ceq):.6e} mol/kg")
    print(f"   fallback faces (Ceq > 1e-2, <= 0 or non-finite): {len(bad)}")

    cache = {}
    def ref(T):
        k = round(T, 6)
        if k not in cache:
            r = vc.evaluate("x", k, mA)
            cache[k] = (float(r["Ceq"]), float(r["Ceq_if_k2_div_a"]))
        return cache[k]

    good = [i for i in range(len(ceq)) if i not in set(bad)]
    if Tl is not None:
        rel = [(ceq[i] - ref(Tl[i])[0]) / ref(Tl[i])[0] for i in good]
        rel_div = [(ceq[i] - ref(Tl[i])[1]) / ref(Tl[i])[1] for i in good]
        worst = max(abs(x) for x in rel)
        n_out = sum(abs(x) > tol for x in rel)
        print(f"   face-by-face vs independent solution: max |rel. diff| = {worst:.2e}, "
              f"mean = {sum(rel)/len(rel):+.2e}, faces outside tolerance = {n_out}")
        print(f"   (same comparison against the k2/a form: max |rel. diff| = "
              f"{max(abs(x) for x in rel_div):.2e})")
    else:
        grid = [Trange[0] + (Trange[1] - Trange[0]) * k / 20 for k in range(21)]
        lo = min(ref(t)[0] for t in grid) * (1 - tol)
        hi = max(ref(t)[0] for t in grid) * (1 + tol)
        n_out = sum(not (lo <= ceq[i] <= hi) for i in good)
        mid = ref(sum(Trange) / 2)[0]
        print(f"   independent Ceq over that T range: {lo:.6e} - {hi:.6e} mol/kg "
              f"(incl. tolerance); faces outside = {n_out}")
        print(f"   mean CeqFace / reference at mid-range T = "
              f"{sum(ceq[i] for i in good)/len(good)/mid:.5f}")
        print("   (add the mesh owner and boundary files for an exact face-by-face check)")
    ok = n_out == 0 and not bad
    print(f"   RESULT: {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=".")
    p.add_argument("--patch", help="override patchName from FACProperties")
    p.add_argument("--tol", type=float, help="relative tolerance (default: set from written precision)")
    args = p.parse_args()

    kw0, ph0 = list(vc.KW_TABLE), list(vc.PH2_TABLE)
    cases = sorted(e for e in os.listdir(args.root)
                   if os.path.isfile(os.path.join(args.root, e, "CeqFace")))
    if not cases:
        sys.exit(f"no <case>/CeqFace folders under {args.root}")
    results = {}
    for c in cases:
        vc.KW_TABLE[:], vc.PH2_TABLE[:] = kw0, ph0
        try:
            results[c] = check_case(os.path.join(args.root, c), args.patch, args.tol)
        except Exception as e:
            print(f"\n== {c}: ERROR: {e}")
            results[c] = False
    print("\nSummary: " + ", ".join(f"{c} {'PASS' if ok else 'FAIL'}" for c, ok in results.items()))
    sys.exit(0 if all(results.values()) else 1)


if __name__ == "__main__":
    main()
