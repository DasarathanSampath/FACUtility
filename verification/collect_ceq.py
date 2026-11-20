#!/usr/bin/env python3
"""
Collect CeqFace values from FACUtility verification runs into the CSV that
verify_chemistry.py --compare reads.

Usage
  python3 collect_ceq.py --map runs.txt --patch downstream --out facutility_results.csv

runs.txt has one line per run:   <label> <case directory> [time directory]
  e.g.  T423.15_case12_pH2table   ~/verif/T423_case12      0
        T473.15_case12_a0.01      ~/verif/T473_a001
If the time directory is omitted, the latest numeric time directory is used.
Lines starting with # are ignored.

For each run, the script reads <case>/<time>/CeqFace (the surfaceScalarField
FACUtility writes), takes the values on the selected patch, checks that they
are uniform (the test runs use uniform T), and writes their mean. Decomposed
cases are read from processor*/<time>/CeqFace if the reconstructed file is
missing.
"""

import argparse
import glob
import os
import re
import statistics
import sys


def latest_time(d):
    times = []
    for name in os.listdir(d):
        try:
            v = float(name)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(d, name)):
            times.append((v, name))
    if not times:
        raise FileNotFoundError(f"no time directories in {d}")
    return max(times)[1]


def strip_comments(text):
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def patch_values(path, patch):
    """Return (values, raw_tokens) on `patch` from an OpenFOAM field file."""
    text = strip_comments(open(path).read())
    bf = text.find("boundaryField")
    if bf < 0:
        raise ValueError(f"{path}: no boundaryField")
    m = re.search(r"\b" + re.escape(patch) + r"\s*\{", text[bf:])
    if not m:
        raise ValueError(f"{path}: patch '{patch}' not found")
    start = bf + m.end()
    depth, i = 1, start
    while depth:
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        i += 1
    block = text[start:i - 1]
    vm = re.search(r"\bvalue\s+(uniform|nonuniform)\s+(.*?);", block, flags=re.S)
    if not vm:
        return [], []  # empty patch on this processor
    kind, body = vm.group(1), vm.group(2).strip()
    if kind == "uniform":
        return [float(body)], [body]
    lm = re.search(r"List<scalar>\s*(\d+)?\s*\((.*)\)", body, flags=re.S)
    if not lm:
        lm = re.search(r"List<scalar>\s*(\d+)\s*\{\s*(\S+)\s*\}", body)
        if lm:  # n{value} shorthand
            return [float(lm.group(2))] * int(lm.group(1)), [lm.group(2)]
        raise ValueError(f"{path}: cannot parse value on patch '{patch}'")
    toks = lm.group(2).split()
    return [float(t) for t in toks], toks


def sig_digits(tok):
    mant = re.split(r"[eE]", tok)[0].lstrip("+-").replace(".", "").lstrip("0")
    return len(mant)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--map", help="text file: label  case_dir|CeqFace_file  [time]")
    src.add_argument("--dir", help="folder of copied CeqFace files, one per label (see README)")
    p.add_argument("--patch", required=True, help="patch name, as patchName in FACProperties")
    p.add_argument("--field", default="CeqFace")
    p.add_argument("--out", default="facutility_results.csv")
    args = p.parse_args()

    rows, problems = [], 0
    for parts in _entries(args):
        label, case = parts[0], os.path.expanduser(parts[1])
        try:
            files = []
            t = parts[2] if len(parts) > 2 else None
            if os.path.isfile(case):
                files = [case]
            elif glob.glob(os.path.join(case, "processor*", args.field)):
                files = sorted(glob.glob(os.path.join(case, "processor*", args.field)))
            elif os.path.isfile(os.path.join(case, args.field)):
                files = [os.path.join(case, args.field)]
            elif os.path.isdir(case) and (t or any(_isnum(n) for n in os.listdir(case))):
                tt = t or latest_time(case)
                f = os.path.join(case, tt, args.field)
                if os.path.isfile(f):
                    files = [f]
            if not files:
                procs = sorted(glob.glob(os.path.join(case, "processor*")))
                for pd in procs:
                    tt = t or latest_time(pd)
                    f = os.path.join(pd, tt, args.field)
                    if os.path.isfile(f):
                        files.append(f)
            if not files:
                raise FileNotFoundError(f"{args.field} not found under {case}")

            vals, toks = [], []
            for f in files:
                v, tk = patch_values(f, args.patch)
                vals += v
                toks += tk
            if not vals:
                raise ValueError(f"no values on patch '{args.patch}'")

            mean = statistics.fmean(vals)
            spread = (max(vals) - min(vals)) / abs(mean) if mean else 0.0
            digits = max(sig_digits(x) for x in toks)
            note = []
            if spread > 1e-9:
                note.append(f"NOT UNIFORM (spread {spread:.2e}) - check T is uniform")
                problems += 1
            if digits < 8:
                note.append(f"only {digits} significant digits - raise writePrecision")
                problems += 1
            print(f"  {label:30s} n={len(vals):6d}  Ceq={mean:.10e}  {'; '.join(note) or 'ok'}")
            rows.append((label, mean))
        except Exception as e:
            print(f"  {label:30s} ERROR: {e}")
            problems += 1

    with open(args.out, "w") as f:
        f.write("label,Ceq_FACUtility\n")
        for label, v in rows:
            f.write(f"{label},{v:.12e}\n")
    print(f"\n{len(rows)} runs written to {args.out}; {problems} warning(s)")
    sys.exit(1 if problems else 0)


def _entries(args):
    """Yield [label, path, (time)] from --map or --dir."""
    if args.map:
        for line in open(args.map):
            line = line.split("#", 1)[0].strip()
            if line:
                yield line.split()
        return
    d = os.path.expanduser(args.dir)
    for name in sorted(os.listdir(d)):
        path = os.path.join(d, name)
        if name.startswith(".") or name.endswith((".csv", ".py", ".txt", ".md")):
            continue
        label = name[:-len("." + args.field)] if name.endswith("." + args.field) else name
        yield [label, path]


def _isnum(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main()
