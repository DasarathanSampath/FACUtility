# Files added to the FACUtility repository for this release

| Folder | What it holds | Paper |
| --- | --- | --- |
| `verification/` | Independent check of the chemistry solver | Section 4.7.1, Tables 2-3 |
| `postProcessingScripts/` | Scripts that make the validation tables | Tables 4-11 |
| `comparisonData/measured/` | Measured values from the published papers | Tables 6-11 |
| `comparisonData/simulated/` | Values extracted from the FACUtility runs | Tables 3-11, Section 4.6.1 |

Requirements: Python 3.8 or newer, with `numpy`, `pandas` and `mpmath`.

    pip install numpy pandas mpmath

## verification/

See `verification/README.md`.

- `verify_chemistry.py`: solves Eqs. (6)-(9) on its own, with no FACUtility code, and writes `reference_values.csv` (21 conditions). Table 2 is taken from this file.
- `verify_cases.py`: checks every face of `CeqFace` from the validation runs against the independent solution, and counts failed solves. Table 3 is made with this script.
- `collect_ceq.py`: collects `CeqFace` values from FACUtility runs into one CSV.
- `selftest_correct.csv` and `selftest_k2div.csv`: quick self-test of the comparison. The first must pass and the second must fail:

      python3 verify_chemistry.py --compare selftest_correct.csv
      python3 verify_chemistry.py --compare selftest_k2div.csv

## postProcessingScripts/

| Script | Input | Output | Paper |
| --- | --- | --- | --- |
| `mesh_window_stats.py` | ParaView CSV of the wall patch, one per mesh | TL and y+ statistics over the 5-20 mm window | Tables 4-7 |
| `process_elbow_slices_4.py` | ParaView CSV of the five slices | TL at the 40 sample points | Tables 8-9, Fig. 12 |
| `case3_compare.py` | output of `process_elbow_slices_4.py` | normalized P, M, E and the statistics | Tables 8-10 |
| `case4_sensors.py` | ParaView CSV of the elbow wall | sensor-area averages and the calibrated chromium factor | Table 11 |

The ParaView export steps are given at the top of each script.

## comparisonData/measured/

| File | Source | Paper |
| --- | --- | --- |
| `case1_case2_yoneda2016_run07.csv` | Yoneda et al. (2016), run #07, Figs. 8(c) and 13(a) | Tables 6-7 |
| `case3_yoneda2018_normalized.csv` | Yoneda et al. (2018), Table 7 / Fig. 16 | Tables 8-9 |
| `case4_utanohara2019_run3.csv` | Utanohara and Murase (2019), Run 3 | Table 11 |

## comparisonData/simulated/

See `comparisonData/simulated/README.md`.
