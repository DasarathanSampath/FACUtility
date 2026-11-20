#!/usr/bin/env python3
"""
Cases 1-2: thinning (TL) and y+ statistics over the 5-20 mm comparison window
for each mesh or inflation setting (Tables 4-7 of the paper).

One CSV per mesh, exported from ParaView from the wall patch used for the
comparison (e.g. 'downstream'):
  reader: tick only that patch > Filters > Cell Centers (Vertex Cells on)
  > File > Save Data > CSV  (keep Points:0..2, TLPatch, yPlus)

Usage
  python3 mesh_window_stats.py A.csv B.csv C.csv D.csv E.csv F.csv
  python3 mesh_window_stats.py inf05.csv inf08.csv --span-coord Points:2

Options
  --coord        streamwise coordinate column (default Points:0)
  --origin       coordinate of the window origin, i.e. the upstream edge of
                 the specimen [m]; default: minimum of --coord on the patch
  --window       start and end of the window [mm] (default 5 20)
  --band         target y+ band (default 30 50)
  --span-coord   3D cases: spanwise coordinate; adds a centreline-only mean
                 (faces within --centre-width of the span midpoint)
  --centre-width full width of the centreline strip [mm] (default 2)

Output: one row per file with
  TL mean / max / min in the window, the location of the max,
  y+ min / 10th / mean / 90th / max in the window,
  % of window faces with y+ inside the band,
  and (3D) the centreline-only TL mean.
Written to mesh_window_stats.csv.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path):
    """ParaView CSV, or FACUtility .dat (whitespace, '# columns:' header)."""
    with open(path) as fh:
        head = [next(fh, "") for _ in range(20)]
    if any(h.lower().startswith("# columns:") for h in head):
        cols = next(h for h in head if h.lower().startswith("# columns:")).split(":", 1)[1].split()
        d = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=cols)
        t = [h.split("=", 1)[1].strip() for h in head if h.startswith("# time")]
        d.attrs["time"] = t[0] if t else ""
        return d
    return pd.read_csv(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+")
    ap.add_argument("--coord", default="Points:0")
    ap.add_argument("--origin", type=float)
    ap.add_argument("--window", type=float, nargs=2, default=[5.0, 20.0])
    ap.add_argument("--band", type=float, nargs=2, default=[30.0, 50.0])
    ap.add_argument("--span-coord")
    ap.add_argument("--centre-width", type=float, default=2.0)
    ap.add_argument("--out", default="mesh_window_stats.csv")
    ap.add_argument("--tl-col", help="TL column (default: TLPatch, else TL)")
    a = ap.parse_args()

    rows = []
    for f in a.csv:
        d = read_table(f)
        tlc = a.tl_col or ("TLPatch" if "TLPatch" in d else "TL")
        if a.coord not in d and a.coord == "Points:0" and "x" in d:
            a.coord = "x"
        if a.span_coord == "Points:2" and "Points:2" not in d and "z" in d:
            a.span_coord = "z"
        x = d[a.coord].to_numpy(float)
        x0 = a.origin if a.origin is not None else x.min()
        xmm = (x - x0) * 1000.0
        w = (xmm >= a.window[0]) & (xmm <= a.window[1])
        tl = d[tlc].to_numpy(float)[w]
        yp = d["yPlus"].to_numpy(float)[w]
        if not w.any():
            raise SystemExit(f"{f}: no faces in the window - check --coord / --origin")
        imax = np.argmax(tl)
        r = {
            "mesh": Path(f).stem, "time": d.attrs.get("time", ""), "origin_m": x0,
            "n_faces": int(w.sum()),
            "TL_mean": tl.mean(), "TL_max": tl.max(),
            "TL_max_at_mm": xmm[w][imax], "TL_min": tl.min(),
            "yplus_min": yp.min(), "yplus_p10": np.percentile(yp, 10),
            "yplus_mean": yp.mean(), "yplus_p90": np.percentile(yp, 90),
            "yplus_max": yp.max(),
            "pct_in_band": 100.0 * np.mean((yp >= a.band[0]) & (yp <= a.band[1])),
        }
        if a.span_coord:
            z = d[a.span_coord].to_numpy(float)[w]
            zc = 0.5 * (z.min() + z.max())
            c = np.abs(z - zc) * 1000.0 <= a.centre_width / 2
            r["TL_mean_centreline"] = tl[c].mean() if c.any() else np.nan
            r["n_centreline"] = int(c.sum())
        rows.append(r)

    t = pd.DataFrame(rows)
    t.to_csv(a.out, index=False)
    pd.set_option("display.width", 200)
    print(f"Window {a.window[0]:g}-{a.window[1]:g} mm from the upstream edge; "
          f"y+ band {a.band[0]:g}-{a.band[1]:g}")
    print(t.round(4).to_string(index=False))
    if len(t) > 1:
        print(f"\nTL mean across these meshes: {t.TL_mean.min():.4f} - {t.TL_mean.max():.4f} "
              f"(spread {100*(t.TL_mean.max()/t.TL_mean.min()-1):.0f}%)")
    print(f"Written: {a.out}")


if __name__ == "__main__":
    main()
