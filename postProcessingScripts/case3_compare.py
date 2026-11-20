#!/usr/bin/env python3
"""
Case 3 (Yoneda et al. 2018 elbow): build Tables 8-10 from the output of
process_elbow_slices_4.py.

Inputs (from process_elbow_slices_4.py, folder processed_slices/):
  TLPatch_by_plane_and_angle.csv   predicted values at 0, 45, ..., 315 deg
                                   (rows Plane A..E)
  *_with_angles.csv                all wall points per slice (for the
                                   predicted ring maximum)

Normalisation: the predicted values are divided by the
maximum of the same 40 sampled points, matching the measured values, which
Yoneda et al. normalise by the maximum of their 40 measured points. The
largest predicted entry is therefore exactly 1. The maximum over all wall
points on the five slices is reported separately, with its location.

Outputs (processed_slices/):
  case3_comparison_unrounded.csv   P, M, E for all 40 points (unrounded)
  case3_tables_8_9.csv             P, M, E rounded as in Tables 8-9
  case3_statistics.csv             Table 10: aggregate and per-plane metrics
  plus a printed summary of the measured and predicted peak locations.

Conventions
  E    = (M - P) / M      (as in Tables 6-7 of the submitted manuscript)
  bias = mean(P - M)      (Table 10)
  r2   = squared Pearson correlation of P and M
  R2   = coefficient of determination, 1 - SS_res / SS_tot (also reported)
"""

from pathlib import Path

import numpy as np
import pandas as pd

FOLDER = Path("processed_slices")
ANGLES = [0, 45, 90, 135, 180, 225, 270, 315]
PLANES = ["A", "B", "C", "D", "E"]

# Measured values, normalised by the measured maximum, from Yoneda et al.
# (2018), Table 7 / Fig. 16, as in Tables 8-9 of the manuscript.
MEASURED = {
    "A": [0.48, 0.30, 0.40, 0.92, 1.00, 1.00, 0.52, 0.38],
    "B": [0.24, 0.20, 0.38, 0.62, 0.82, 0.84, 0.40, 0.24],
    "C": [0.08, 0.17, 0.38, 0.64, 0.54, 0.60, 0.38, 0.20],
    "D": [0.34, 0.32, 0.40, 0.70, 0.24, 0.54, 0.50, 0.32],
    "E": [0.50, 0.50, 0.79, 0.70, 0.18, 0.82, 0.50, 0.48],
}


def plane_key(label):
    return str(label).replace("Plane", "").strip().upper()


def main():
    wide = pd.read_csv(FOLDER / "TLPatch_by_plane_and_angle.csv", index_col=0)
    wide.index = [plane_key(i) for i in wide.index]
    wide.columns = [int(str(c).replace("°", "").strip()) for c in wide.columns]
    wide = wide.loc[PLANES, ANGLES]

    sampled_max = wide.to_numpy().max()
    P = wide / sampled_max
    M = pd.DataFrame(MEASURED, index=ANGLES).T.loc[PLANES]

    # Long, unrounded comparison table
    rows = []
    for pl in PLANES:
        for a in ANGLES:
            p, m = P.loc[pl, a], M.loc[pl, a]
            rows.append({"Plane": pl, "Angle_deg": a,
                         "P_sampled_norm": p,
                         "P_as_exported": wide.loc[pl, a],
                         "M": m, "E": (m - p) / m, "P_minus_M": p - m})
    long = pd.DataFrame(rows)
    long.to_csv(FOLDER / "case3_comparison_unrounded.csv", index=False)

    # Tables 6-7 layout, rounded
    t = []
    for pl in PLANES:
        r = {"Plane": pl}
        for a in ANGLES:
            p, m = P.loc[pl, a], M.loc[pl, a]
            r[f"{a}_P"] = round(p, 2)
            r[f"{a}_M"] = m
            r[f"{a}_E%"] = round(100 * (m - p) / m)
        t.append(r)
    pd.DataFrame(t).to_csv(FOLDER / "case3_tables_8_9.csv", index=False)

    # Table 10
    def metrics(p, m):
        d = p - m
        return {
            "N": len(d),
            "RMSE": np.sqrt(np.mean(d**2)),
            "MAE": np.mean(np.abs(d)),
            "Bias_P_minus_M": np.mean(d),
            "MAPE_%": 100 * np.mean(np.abs(d) / m),
            "NRMSE_%": 100 * np.sqrt(np.mean(d**2)) / np.mean(m),
            "r2_squared_Pearson": np.corrcoef(p, m)[0, 1] ** 2,
            "R2_coeff_determination": 1 - np.sum(d**2) / np.sum((m - m.mean())**2),
        }

    stats = [{"Set": "All planes", **metrics(long.P_sampled_norm.to_numpy(), long.M.to_numpy())}]
    for pl in PLANES:
        s = long[long.Plane == pl]
        stats.append({"Set": f"Plane {pl}", **metrics(s.P_sampled_norm.to_numpy(), s.M.to_numpy())})
    stats = pd.DataFrame(stats)
    stats.to_csv(FOLDER / "case3_statistics.csv", index=False)

    # Peak locations
    mmax = M.to_numpy().max()
    meas_peaks = [(pl, a) for pl in PLANES for a in ANGLES if M.loc[pl, a] == mmax]
    pi = np.unravel_index(np.argmax(P.to_numpy()), P.shape)
    pred_peak = (PLANES[pi[0]], ANGLES[pi[1]])

    ring_max, ring_loc = -np.inf, None
    for f in sorted(FOLDER.glob("*_with_angles.csv")):
        d = pd.read_csv(f)
        # TLPatch_by_plane_and_angle.csv holds values already divided by the
        # slice-wide maximum, so compare in the same units.
        col = "TLPatch_Normalized" if "TLPatch_Normalized" in d.columns else "TLPatch"
        i = d[col].idxmax()
        if d.loc[i, col] > ring_max:
            ring_max = d.loc[i, col]
            ring_loc = (f.stem.replace("_with_angles", ""), d.loc[i, "Circumferential_Angle_deg"])

    pd.set_option("display.width", 160)
    print("Measured peak (M = 1):", ", ".join(f"Plane {p} {a} deg" for p, a in meas_peaks))
    print(f"Predicted peak among the 40 sampled points: Plane {pred_peak[0]} {pred_peak[1]} deg")
    if ring_loc:
        print(f"Predicted maximum over all slice wall points: {ring_loc[0]} at "
              f"{ring_loc[1]:.1f} deg = {ring_max / sampled_max:.4f} x the sampled maximum")
    print(f"Sampled maximum TLPatch (normalising value): {sampled_max:.6g}")
    print("\nPredicted at the measured peak(s):",
          ", ".join(f"Plane {p} {a} deg: P = {P.loc[p, a]:.3f}" for p, a in meas_peaks))
    print("\nTable 10:")
    print(stats.round(3).to_string(index=False))
    print(f"\nWritten to {FOLDER.resolve()}")


if __name__ == "__main__":
    main()
