#!/usr/bin/env python3
"""
Case 1 (Mesh F): face-by-face comparison of the thinning rate from the
gradient-based and the traction-based shear velocity (Section 4.6.1, Eq. 15).

For each face, Eq. (15) gives
    TL_t = TL_g * (Utau_t / Utau_g)^2
because all other terms of km are the same for both methods.

Input CSV columns: Xdistance [m], Gradient_Utau, Traction_Utau, Gradient_TL.

Usage
  python3 case1_traction_vs_gradient.py case1_meshF_TL_traction_vs_gradient.csv
  (--origin: x of the upstream edge of the specimen, default 0.253 m;
   --window: start and end of the comparison window in mm, default 5 20)

All means are arithmetic means over the faces. The faces of the 2D patch
have the same size, so this equals an area-weighted mean.
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--origin", type=float, default=0.253)
    ap.add_argument("--window", type=float, nargs=2, default=[5.0, 20.0])
    ap.add_argument("--out", default="case1_meshF_TL_paired_window.csv")
    a = ap.parse_args()

    d = pd.read_csv(a.csv, encoding="utf-8-sig")
    d.columns = [c.strip() for c in d.columns]
    s_mm = (d["Xdistance"] - a.origin) * 1000.0
    w = d[(s_mm >= a.window[0] - 1e-9) & (s_mm <= a.window[1] + 1e-9)].copy()
    w["s_mm"] = s_mm[w.index]
    w["TL_g"] = w["Gradient_TL"]
    w["TL_t"] = w["TL_g"] * (w["Traction_Utau"] / w["Gradient_Utau"]) ** 2
    w["diff"] = w["TL_t"] - w["TL_g"]
    w["rel_diff_pct"] = 100.0 * w["diff"] / w["TL_g"]

    n = len(w)
    imax = w["diff"].abs().idxmax()
    print(f"Window {a.window[0]:g}-{a.window[1]:g} mm (origin x = {a.origin} m): {n} faces")
    print(f"  mean of TL_t - TL_g   = {w['diff'].mean():.2e} mm/10,000 h")
    print(f"  RMSE                  = {np.sqrt((w['diff']**2).mean()):.2e} mm/10,000 h")
    print(f"  MAE                   = {w['diff'].abs().mean():.2e} mm/10,000 h")
    print(f"  largest |difference|  = {abs(w.loc[imax, 'diff']):.2e} at {w.loc[imax, 's_mm']:.1f} mm"
          f" (TL_g there = {w.loc[imax, 'TL_g']:.3f})")
    print(f"  faces with |rel. difference| > 5% = {(w['rel_diff_pct'].abs() > 5).sum()}")
    print(f"  median |rel. difference|          = {w['rel_diff_pct'].abs().median():.3f} %")
    mg, mt = w["TL_g"].mean(), w["TL_t"].mean()
    print(f"  window-mean TL: gradient {mg:.4f} -> traction {mt:.4f} ({100*(mt/mg-1):+.3f} %)")
    w[["Xdistance", "s_mm", "Gradient_Utau", "Traction_Utau", "TL_g", "TL_t", "diff",
       "rel_diff_pct"]].to_csv(a.out, index=False)
    print(f"Written: {a.out}")


if __name__ == "__main__":
    main()
