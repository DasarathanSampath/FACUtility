# FACUtility chemistry verification test

Independent check of `FACChemistryModel` against Eqs. (6)-(9), with no shared code.
Requires Python 3.8+ and `mpmath` (`pip install mpmath`).

Files
- `verify_chemistry.py`  independent solver; writes reference values and compares
- `collect_ceq.py`       reads CeqFace from FACUtility runs into a CSV
- `reference_values.csv` reference output for 21 conditions

## 1. Reference values

    python3 verify_chemistry.py

Writes `reference_values.csv`: for 21 conditions (T = 25-200 °C, three alkali levels,
and p(H2) overridden to 1e-3 and 1e-6 atm so that a = 0.1 and 0.01) it lists K0-K3,
m(H+), pH at temperature, the four species molalities and the total Ceq. The
hydrogen-ion balance is solved twice at 50-digit precision (cubic roots and bisection),
and the two solutions must agree to 1e-30. `Ceq_if_k2_div_a` gives the value a code
would produce if Fe(OH)2 were divided by a, so a failure can be diagnosed.

## 2. FACUtility runs

For each row to test, copy any validation case (CFD fields already converged) and set:

- `system/controlDict`: `writePrecision 12;`
- the time directory FACUtility reads: `T` uniform at `T_K`, on the internal field
  and on the evaluated patch (`internalField uniform 423.15;` and
  `value uniform 423.15;` with a fixedValue/calculated patch type)
- `system/FACProperties`: `mA` = the row's `mA`
- a = 0.1 / 0.01 rows only: `constant/pH2_vs_T.dat` replaced by a two-row table
  with `pH2_atm` at 273.15 K and 673.15 K

Run FACUtility. `CeqFace` in the time directory holds the total molality (mol/kg)
from Listing 8, line 48, before the Eq. 6 volume conversion done in FACThinningModel.

## 3. Collect

**Option A: copied CeqFace files only (no full cases needed).** Put the scripts,
`reference_values.csv` and the copied files in one folder, named by label:

    Ceq_Verification/
      T423.15_case12_pH2table/CeqFace             # a subfolder per label, or
      T473.15_case12_a0.01.CeqFace                # a file named <label>.CeqFace
      T423.15_case12_a0.1/processor0/CeqFace      # decomposed: one per processor
      T423.15_case12_a0.1/processor1/CeqFace

then, from inside that folder

    python3 collect_ceq.py --dir . --patch downstream

**Option B: full case directories.** List the runs in a text file, one per line: `label  case_directory  [time]`

    T423.15_case12_pH2table   /path/to/T423_case12
    T473.15_case12_a0.01      /path/to/T473_a001   0

then

    python3 collect_ceq.py --map runs.txt --patch downstream --out facutility_results.csv

It warns if the values on the patch are not uniform (T was not uniform) or have fewer
than 8 significant digits (writePrecision too low).

## 4. Compare

    python3 verify_chemistry.py --compare facutility_results.csv

Pass criterion: relative difference <= 1e-6 (`--tol` to change). Exit status 1 on
any failure, so it can run in CI. Case4_149.2C is the only temperature between
table nodes; a small failure there alone points to the table interpolation scheme
(linear in log10 Kw is assumed), not the chemistry.

## Checking the validation runs themselves (Cases 1-4) and nFail

`verify_cases.py` compares every face of CeqFace from the validation runs against
the independent solution at that face's temperature, and counts faces holding the
Newton-Raphson fallback value.

    nFail_verification/
      Kw_vs_T.dat, pH2_vs_T.dat   constant/ tables (here, or inside a case folder)
      case1/CeqFace               <time>/CeqFace
      case1/T                     <time>/T  (same time directory)
      case1/FACProperties         system/FACProperties  (mA and patchName read from it)
      case1/owner                 optional: constant/polyMesh/owner (ascii)
      case1/boundary              optional: constant/polyMesh/boundary

Run `python3 verify_cases.py` from that folder.

FACUtility evaluates each face at its owner-cell temperature T[own] (Listing 8,
line 22), so the wall patch type in T does not matter. With owner and boundary the
check is exact, face by face. Without them, each face is checked against the range
of Ceq spanned by the lowest and highest cell temperature in the domain; that bound
is tight when T is nearly uniform. A failed solve clamps m(H+) to 1e-20, giving
Ceq ~ 1e4 mol/kg, so fallback faces cannot be missed. The tolerance is set from the
written precision (5e-6 for writePrecision 6).
