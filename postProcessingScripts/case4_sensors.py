#!/usr/bin/env python3
"""
Case 4 (Utanohara & Murase 2019, Run 3): sensor-footprint averages of TLPatch
from a ParaView export of the elbow wall surface, and the chromium factor that
makes the mean sim/exp over the eight straight-section sensors equal to 1.

ParaView export (cell data, one row per wall face)
  1. Reader: tick only the wall patches (internalMesh off). Apply.
  2. Filters > Cell Size        (adds 'Area'; tick only "Compute Area")
  3. Filters > Cell Centers     (tick "Vertex Cells"; keeps cell data)
  4. Select CellCenters1 > File > Save Data > CSV
     Keep Points:0..2, TLPatch, Area, yPlus.

Usage
  python3 case4_sensors.py walls.csv --factor 0.112
  (--factor = the chromium factor of the run the CSV comes from)

Geometry (defaults from the CAD model, same as the test section):
  elbow entrance at x0 = 0.0064 m on the upstream centreline y = yc = 0.1223 m,
  two 45 deg bends of centreline radius R = 0.0762 m joined by a straight
  section L = 0.054 m, D = 0.0495 m; flow turns from +x to -y in the z = 0
  plane. Positions are located by centreline distance s from the entrance:
      upstream     s = -9.35 D
      elbow        mid-length of the straight section between the bends
      downstream   2.68 D and 3.68 D after the elbow exit
  (x/D = -9.35, 1.75, 6.2, 7.2 in the paper's Fig. 7 convention).
  Circumferential positions: extrados = outer side of the bend (+y upstream,
  +x downstream), intrados = opposite, upper = +z, lower = -z.

Footprint: 20 mm along the centreline x 3 mm around the circumference, as the
sensor plates of Utanohara and Murase (2019). Area-weighted mean of TLPatch. Faces
with yPlus > 70 are reported separately: each average is given
with and without them.
"""

import argparse

import numpy as np
import pandas as pd

D = 0.0495
EXP = {  # Run 3 measured FAC rate, from Table 11 of the manuscript
    ("upstream", "extrados"): 0.134, ("upstream", "intrados"): 0.110,
    ("elbow", "extrados"): 0.353, ("elbow", "intrados"): 0.271, ("elbow", "upper"): 0.298,
    ("down6.2", "intrados"): 0.165, ("down6.2", "upper"): 0.170, ("down6.2", "lower"): 0.200,
    ("down7.2", "extrados"): 0.172, ("down7.2", "intrados"): 0.100, ("down7.2", "upper"): 0.207,
}
XD = {"upstream": -9.35, "elbow": 1.75, "down6.2": 6.2, "down7.2": 7.2}


def centreline(x0, yc, R, L, n=4000):
    """Dense centreline points, arc length s (0 at entrance), unit tangents."""
    c = np.cos(np.pi / 4)
    up = np.c_[np.linspace(x0 - 0.6, x0, n), np.full(n, yc)]
    t = np.linspace(0, np.pi / 4, n)
    a1 = np.c_[x0 + R * np.sin(t), yc - R * (1 - np.cos(t))]
    st = a1[-1] + np.linspace(0, L, n)[:, None] * np.array([c, -c])
    C2 = st[-1] + R * np.array([-c, -c])
    a0 = np.arctan2(st[-1, 1] - C2[1], st[-1, 0] - C2[0])
    t2 = np.linspace(a0, a0 - np.pi / 4, n)
    a2 = np.c_[C2[0] + R * np.cos(t2), C2[1] + R * np.sin(t2)]
    dn = np.c_[np.full(n, a2[-1, 0]), np.linspace(a2[-1, 1], a2[-1, 1] - 0.6, n)]
    P = np.vstack([up, a1[1:], st[1:], a2[1:], dn[1:]])
    seg = np.diff(P, axis=0)
    s = np.r_[0, np.cumsum(np.hypot(*seg.T))]
    s -= s[n - 1]                                    # entrance at s = 0
    tan = np.vstack([seg, seg[-1]])
    tan /= np.linalg.norm(tan, axis=1)[:, None]
    s_exit = 2 * R * np.pi / 4 + L
    return P, s, tan, s_exit


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--factor", type=float, required=True,
                    help="chromium factor of the run the CSV comes from")
    ap.add_argument("--x0", type=float, default=0.0064)
    ap.add_argument("--yc", type=float, default=0.1223)
    ap.add_argument("--zc", type=float, default=0.0)
    ap.add_argument("--R", type=float, default=0.0762)
    ap.add_argument("--L", type=float, default=0.054)
    ap.add_argument("--axial", type=float, default=0.020, help="footprint length [m]")
    ap.add_argument("--width", type=float, default=0.003, help="footprint width [m]")
    ap.add_argument("--yplus-max", type=float, default=70.0)
    ap.add_argument("--out", default="case4_sensor_values.csv")
    a = ap.parse_args()

    df = pd.read_csv(a.csv)
    X = df[["Points:0", "Points:1", "Points:2"]].to_numpy(float)
    tl = df["TLPatch"].to_numpy(float)
    area = df["Area"].to_numpy(float) if "Area" in df else np.ones(len(df))
    yp = df["yPlus"].to_numpy(float) if "yPlus" in df else np.zeros(len(df))
    if "Area" not in df:
        print("note: no 'Area' column - using unweighted means")

    P, s_cl, tan, s_exit = centreline(a.x0, a.yc, a.R, a.L)
    from scipy.spatial import cKDTree
    _, idx = cKDTree(P).query(X[:, :2])
    s = s_cl[idx]
    rel = np.c_[X[:, :2] - P[idx], X[:, 2] - a.zc]
    # remove the small axial component so rel is a pure radial vector
    t3 = np.c_[tan[idx], np.zeros(len(idx))]
    rel -= (np.sum(rel * t3, axis=1))[:, None] * t3
    out_dir = np.c_[-tan[idx, 1], tan[idx, 0], np.zeros(len(idx))]   # extrados (left normal)
    phi = np.degrees(np.arctan2(rel[:, 2], np.sum(rel * out_dir, axis=1))) % 360  # 0=extrados, 90=upper(+z)
    rad = np.linalg.norm(rel, axis=1)

    s_station = {"upstream": -9.35 * D, "elbow": a.R * np.pi / 4 + a.L / 2,
                 "down6.2": s_exit + 2.68 * D, "down7.2": s_exit + 3.68 * D}
    phi_pos = {"extrados": 0.0, "upper": 90.0, "intrados": 180.0, "lower": 270.0}

    rows = []
    for (st, pos), exp in EXP.items():
        dphi = (phi - phi_pos[pos] + 180) % 360 - 180
        arc = np.radians(dphi) * (D / 2)
        sel = (np.abs(s - s_station[st]) <= a.axial / 2) & (np.abs(arc) <= a.width / 2)
        n_all = int(sel.sum())
        widened = False
        if n_all < 3:               # footprint narrower than the faces: widen to capture >= 3
            for w in (0.006, 0.010, 0.015):
                sel = (np.abs(s - s_station[st]) <= a.axial / 2) & (np.abs(arc) <= w / 2)
                if sel.sum() >= 3:
                    widened = w
                    break
        keep = sel & (yp <= a.yplus_max)
        mean_all = np.average(tl[sel], weights=area[sel]) if sel.any() else np.nan
        mean_keep = np.average(tl[keep], weights=area[keep]) if keep.any() else np.nan
        rows.append({
            "station": st, "x/D": XD[st], "position": pos, "exp": exp,
            "sim_all_faces": mean_all, "sim_excl_yplus": mean_keep,
            "n_faces": int(sel.sum()), "n_excluded": int((sel & ~keep).sum()),
            "footprint_width_mm": 1000 * (widened or a.width),
            "mean_wall_radius_mm": 1000 * rad[sel].mean() if sel.any() else np.nan,
        })
    r = pd.DataFrame(rows)
    r["sim/exp"] = r["sim_excl_yplus"] / r["exp"]
    r.to_csv(a.out, index=False)

    straight = r[r.station != "elbow"]
    mean_ratio = straight["sim/exp"].mean()
    new_factor = a.factor / mean_ratio
    pd.set_option("display.width", 160)
    print(r.round(4).to_string(index=False))
    print(f"\nRun factor {a.factor}: mean sim/exp over the {len(straight)} straight-section sensors = {mean_ratio:.4f}")
    print(f"Factor giving mean sim/exp = 1: {a.factor} / {mean_ratio:.4f} = {new_factor:.4f}")
    el = r[r.station == "elbow"]["sim/exp"] * new_factor / a.factor
    print(f"Elbow sim/exp at that factor: {', '.join(f'{v:.3f}' for v in el)} "
          f"(1 - sim/exp: {100*(1-el.max()):.0f}% to {100*(1-el.min()):.0f}%)")
    if (r.n_excluded > 0).any():
        print("Faces with yPlus > %.0f fall inside some footprints - see n_excluded." % a.yplus_max)
    else:
        print("No footprint contains faces with yPlus > %.0f." % a.yplus_max)
    hi = yp > a.yplus_max
    if hi.any():
        d = np.min(np.abs(s[hi][:, None] - np.array([s_station["elbow"]])), axis=0)[0]
        print(f"{hi.sum()} faces with yPlus > {a.yplus_max:.0f}; centreline s = "
              f"{1000*s[hi].min():.1f} to {1000*s[hi].max():.1f} mm from the entrance; "
              f"nearest is {1000*d:.1f} mm ({d/D:.2f} D) from the elbow-sensor centre "
              f"(footprint half-length {1000*a.axial/2:.0f} mm)")
    print(f"Written: {a.out}")


if __name__ == "__main__":
    main()
