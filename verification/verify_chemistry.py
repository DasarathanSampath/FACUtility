#!/usr/bin/env python3
"""
Independent verification of the FACUtility chemistry model (Eqs. 6-9).

This script re-implements the equilibrium chemistry directly from the
governing equations in the manuscript, sharing no code with FACUtility:

  Eq. 7   m_b = K_b * m(H+)^(2-b) * a,   a = p(H2)^(1/3)
          b = 0: Fe2+, 1: FeOH+, 2: Fe(OH)2, 3: Fe(OH)3-
  Eq. 8   ln K_b = (-A_b/T + B_b (ln T - 1) + D_b) / R
  Eq. 9   2 K0 a m^3 + (K1 a + 1) m^2 - mA m - (K3 a + Kw) = 0
  Ceq     = sum of the four species molalities (mol/kg); this is the
            quantity FACUtility stores in CeqFace (Listing 8, line 48)

The hydrogen-ion balance (Eq. 9) is a cubic with exactly one positive
root (one sign change in its coefficients), so it is solved here two
independent ways at 50-digit precision:
  (1) all roots of the cubic polynomial (mpmath.polyroots), positive one kept
  (2) bracketed bisection on [0, 1]
and the two must agree to better than 1e-30 (relative).

Usage
  python3 verify_chemistry.py                 # reference table -> reference_values.csv
  python3 verify_chemistry.py --compare F.csv # compare FACUtility output against it

The comparison CSV needs columns: label, Ceq_FACUtility
(label must match a row of reference_values.csv). Pass criterion:
relative difference <= --tol (default 1e-6).
"""

import argparse
import csv
import sys

import mpmath as mp

mp.mp.dps = 50

# ---------------------------------------------------------------- inputs
# Listing 1: constant/chemistryConstants.dat (Sweeton & Baes)
R = mp.mpf("1.987")  # cal/(mol K)
SPECIES = [  # (name, b, Ab, Bb, Db)
    ("Fe2+", 0, mp.mpf("-26876"), mp.mpf("9.81"), mp.mpf("-81.21")),
    ("FeOH+", 1, mp.mpf("-11733"), mp.mpf("3.35"), mp.mpf("-42.72")),
    ("Fe(OH)2", 2, mp.mpf("4615"), mp.mpf("0.0"), mp.mpf("-23.57")),
    ("Fe(OH)3-", 3, mp.mpf("9045"), mp.mpf("0.0"), mp.mpf("-49.37")),
]

# Listing 3: constant/Kw_vs_T.dat  (values are log10(Kw), see Section 3.3)
KW_TABLE = [
    (298.15, -13.997), (323.15, -13.262), (373.15, -12.216),
    (423.15, -11.65), (473.15, -11.3), (523.15, -11.18),
    (533.15, -11.17), (573.15, -11.19),
]

# Listing 4: constant/pH2_vs_T.dat  (atm)
PH2_TABLE = [
    (298.15, 1.0), (323.15, 1.056), (373.15, 0.951), (423.15, 0.712),
    (473.15, 0.478), (523.15, 0.287), (533.15, 0.264), (573.15, 0.159),
]


def interp(table, T):
    """Linear interpolation, clamped at the table ends."""
    T = mp.mpf(T)
    xs = [mp.mpf(str(x)) for x, _ in table]
    ys = [mp.mpf(str(y)) for _, y in table]
    if T <= xs[0]:
        return ys[0]
    if T >= xs[-1]:
        return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= T <= xs[i + 1]:
            w = (T - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + w * (ys[i + 1] - ys[i])
    raise ValueError(T)


def Kb(Ab, Bb, Db, T):
    """Eq. 8 (natural log)."""
    T = mp.mpf(T)
    return mp.e ** ((-Ab / T + Bb * (mp.log(T) - 1) + Db) / R)


# ---------------------------------------------------------------- solve
def solve_mH(K, a, Kw, mA):
    c3 = 2 * K[0] * a
    c2 = K[1] * a + 1
    c1 = -mA
    c0 = -(K[3] * a + Kw)

    # (1) polynomial roots
    roots = mp.polyroots([c3, c2, c1, c0], maxsteps=200, extraprec=200)
    pos = [mp.re(r) for r in roots if abs(mp.im(r)) < mp.mpf("1e-40") and mp.re(r) > 0]
    if len(pos) != 1:
        raise RuntimeError(f"expected one positive root, got {pos}")
    m_poly = pos[0]

    # (2) bisection
    F = lambda m: ((c3 * m + c2) * m + c1) * m + c0
    lo, hi = mp.mpf(0), mp.mpf(1)
    assert F(lo) < 0 < F(hi)
    for _ in range(400):
        mid = (lo + hi) / 2
        if F(mid) > 0:
            hi = mid
        else:
            lo = mid
    m_bis = (lo + hi) / 2

    rel = abs(m_poly - m_bis) / m_poly
    if rel > mp.mpf("1e-30"):
        raise RuntimeError(f"root methods disagree: {rel}")
    return m_poly, F(m_poly)


def evaluate(label, T, mA, pH2=None):
    T = mp.mpf(str(T))
    mA = mp.mpf(str(mA))
    pH2 = interp(PH2_TABLE, T) if pH2 is None else mp.mpf(str(pH2))
    Kw = mp.mpf(10) ** interp(KW_TABLE, T)
    a = mp.cbrt(pH2)
    K = [Kb(Ab, Bb, Db, T) for _, _, Ab, Bb, Db in SPECIES]

    mH, resid = solve_mH(K, a, Kw, mA)

    # Eq. 7, multiplicative in a for every species
    m = [K[b] * mH ** (2 - b) * a for b in range(4)]
    ceq = sum(m)

    # Alternatives the test must be able to reject
    ceq_k2_div = m[0] + m[1] + K[2] / a + m[3]       # only Fe(OH)2 divided
    ceq_all_div = sum(K[b] * mH ** (2 - b) / a for b in range(4))

    # charge balance check (independent of Eq. 9's algebraic form)
    charge = mH + 2 * m[0] + m[1] + (-mA) - Kw / mH - m[3]

    return {
        "label": label,
        "T_K": T, "T_C": T - mp.mpf("273.15"), "pH2_atm": pH2, "a": a,
        "mA": mA, "log10Kw": mp.log10(Kw),
        "K0": K[0], "K1": K[1], "K2": K[2], "K3": K[3],
        "mH": mH, "pH_hot": -mp.log10(mH),
        "m_Fe2+": m[0], "m_FeOH+": m[1], "m_Fe(OH)2": m[2], "m_Fe(OH)3-": m[3],
        "Ceq": ceq,
        "Ceq_if_k2_div_a": ceq_k2_div,
        "Ceq_if_all_div_a": ceq_all_div,
        "ratio_k2div_to_correct": ceq_k2_div / ceq,
        "Eq9_residual": resid,
        "charge_balance_residual": charge,
    }


# ------------------------------------------------------- test conditions
# Temperatures sit on table nodes so the result does not depend on how
# a code interpolates between nodes; Case 4 (149.2 C) is the exception,
# and checks the interpolation as well.
CONDITIONS = []
for T in (298.15, 373.15, 423.15, 473.15):
    for mA, tag in ((-7.94e-6, "case12"), (-1.585e-5, "case3"), (0.0, "neutral")):
        CONDITIONS.append((f"T{T:.2f}_{tag}_pH2table", T, mA, None))
# a far from unity: p(H2) = 1e-3 atm -> a = 0.1; 1e-6 atm -> a = 0.01
for T in (373.15, 423.15, 473.15):
    for pH2, tag in ((1e-3, "a0.1"), (1e-6, "a0.01")):
        CONDITIONS.append((f"T{T:.2f}_case12_{tag}", T, -7.94e-6, pH2))
# Case-specific operating points
CONDITIONS.append(("Case1_2_150C", 423.15, -7.94e-6, None))
CONDITIONS.append(("Case3_150C", 423.15, -1.585e-5, None))
CONDITIONS.append(("Case4_149.2C", 422.35, -1.494e-7, None))

COLUMNS = ["label", "T_K", "T_C", "pH2_atm", "a", "mA", "log10Kw",
           "K0", "K1", "K2", "K3", "mH", "pH_hot",
           "m_Fe2+", "m_FeOH+", "m_Fe(OH)2", "m_Fe(OH)3-", "Ceq",
           "Ceq_if_k2_div_a", "Ceq_if_all_div_a", "ratio_k2div_to_correct",
           "Eq9_residual", "charge_balance_residual"]


def fmt(v):
    return v if isinstance(v, str) else mp.nstr(v, 15)


def build_reference(path):
    rows = [evaluate(*c) for c in CONDITIONS]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow([fmt(r[c]) for c in COLUMNS])
    return rows


def compare(ref_rows, path, tol):
    ref = {r["label"]: r for r in ref_rows}
    n_pass = n_fail = 0
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            lab = row["label"].strip()
            if lab not in ref:
                print(f"  ?  {lab}: no reference row")
                continue
            got = mp.mpf(row["Ceq_FACUtility"])
            exp = ref[lab]["Ceq"]
            rel = abs(got - exp) / exp
            ok = rel <= tol
            n_pass += ok
            n_fail += not ok
            hint = ""
            if not ok:
                for alt in ("Ceq_if_k2_div_a", "Ceq_if_all_div_a"):
                    if abs(got - ref[lab][alt]) / ref[lab][alt] <= tol:
                        hint = f"  <-- matches {alt}"
            print(f"  {'PASS' if ok else 'FAIL'}  {lab}: FACUtility {mp.nstr(got, 10)}"
                  f"  reference {mp.nstr(exp, 10)}  rel.diff {mp.nstr(rel, 3)}{hint}")
    print(f"\n{n_pass} passed, {n_fail} failed (tolerance {tol:g})")
    return n_fail == 0


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="reference_values.csv")
    p.add_argument("--compare", metavar="CSV")
    p.add_argument("--tol", type=float, default=1e-6)
    args = p.parse_args()

    rows = build_reference(args.out)
    print(f"Reference values for {len(rows)} conditions written to {args.out}\n")
    print(f"{'label':32s} {'a':>7s} {'pH_hot':>7s} {'Ceq (mol/kg)':>13s}"
          f" {'if k2/a':>13s} {'ratio':>7s}")
    for r in rows:
        print(f"{r['label']:32s} {mp.nstr(r['a'], 4):>7s} {mp.nstr(r['pH_hot'], 4):>7s}"
              f" {mp.nstr(r['Ceq'], 6):>13s} {mp.nstr(r['Ceq_if_k2_div_a'], 6):>13s}"
              f" {mp.nstr(r['ratio_k2div_to_correct'], 4):>7s}")
    worst = max(abs(r["charge_balance_residual"]) / r["mH"] for r in rows)
    print(f"\nMax relative charge-balance residual: {mp.nstr(worst, 3)}")

    if args.compare:
        print(f"\nComparing {args.compare}:")
        sys.exit(0 if compare(rows, args.compare, args.tol) else 1)


if __name__ == "__main__":
    main()
