#!/usr/bin/env python3
"""
Figure 9 of the paper: Case 2 thinning rate along the streamwise direction.

Plots TL of each wall face (points), the spanwise mean and the centre-line mean
(2 mm wide strip) at each streamwise position, and the means over the 5-20 mm
window (dotted lines). The window means are the values in Table 7.

    python3 postProcessingScripts/case2_profile_plot.py \
        comparisonData/simulated/case2_inflation05.csv --origin 0.253 \
        --out results/case2_TL_profile.png

Requires numpy, pandas and matplotlib.
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--origin", type=float, default=0.253, help="x of the upstream edge [m]")
    ap.add_argument("--window", type=float, nargs=2, default=[5.0, 20.0], help="[mm]")
    ap.add_argument("--centre-width", type=float, default=2.0, help="[mm]")
    ap.add_argument("--xmax", type=float, default=25.5, help="[mm]")
    ap.add_argument("--out", default="case2_TL_profile.png")
    a = ap.parse_args()

    d = pd.read_csv(a.csv)
    x = (d["Points:0"].to_numpy(float) - a.origin) * 1000.0
    z = d["Points:2"].to_numpy(float)
    tl = d["TLPatch"].to_numpy(float)
    zc = 0.5 * (z.min() + z.max())
    cl = np.abs(z - zc) * 1000.0 <= a.centre_width / 2

    w = (x >= a.window[0]) & (x <= a.window[1])
    m_span, m_cl = tl[w].mean(), tl[w & cl].mean()
    i = np.argmax(tl[w])
    print(f"window {a.window[0]:g}-{a.window[1]:g} mm: spanwise mean {m_span:.4f}, "
          f"centre-line mean {m_cl:.4f}, max {tl[w][i]:.4f} at {x[w][i]:.1f} mm")

    k = x <= a.xmax
    g = pd.DataFrame({"x": np.round(x[k], 4), "tl": tl[k], "cl": cl[k]})
    span = g.groupby("x").tl.mean()
    cent = g[g.cl].groupby("x").tl.mean()

    plt.rcParams.update({"font.family": "serif", "font.size": 11, "axes.linewidth": 0.8})
    blue, red, grey = "#1f4fbf", "#c0392b", "#9a9a9a"
    fig, ax = plt.subplots(figsize=(7.0, 3.9), dpi=300)
    ax.axvspan(a.window[0], a.window[1], color="#f0efec", zorder=0, lw=0)
    ax.scatter(x[k], tl[k], s=3, color=grey, alpha=0.5, lw=0, zorder=1,
               label="$T_L$ of each wall face")
    ax.plot(span.index, span.values, color=blue, lw=2, zorder=3, label="Spanwise mean")
    ax.plot(cent.index, cent.values, color=red, lw=2, ls="--", zorder=3,
            label=f"Centre-line mean ({a.centre_width:g} mm strip)")
    for m, c in ((m_span, blue), (m_cl, red)):
        ax.hlines(m, a.window[0], a.window[1], color=c, lw=1.2, ls=":", zorder=2)
        ax.text(a.window[1] + 0.3, m, f"{m:.3f}", va="center", fontsize=10)
    ax.text(0.5 * sum(a.window), 0.02, f"{a.window[0]:g}-{a.window[1]:g} mm window",
            ha="center", fontsize=9, color="#555555")
    ax.set_xlim(0, a.xmax)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Distance (mm)")
    ax.set_ylabel("$T_L$ (mm/10,000 h)")
    ax.grid(color="#dddddd", lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=9, frameon=False, markerscale=3)
    fig.tight_layout()
    fig.savefig(a.out)
    print(f"Written: {a.out}")


if __name__ == "__main__":
    main()
