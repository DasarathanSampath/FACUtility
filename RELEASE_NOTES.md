# FACUtility v1.1.1

No change to the FACUtility code or the validation cases. All results are the same as in v1.1.0.

- Added the owner and boundary files of the FACUtility patch for Cases 1-4
  (comparisonData/simulated/nFail_verification). With these files,
  verify_cases.py checks each face at the temperature of its owner cell (Table 3).
- Added verification/extract_patch_owner.py, which makes these files from the mesh.
- Added postProcessingScripts/case2_profile_plot.py, which makes Figure 9
  (Case 2: face values, spanwise mean and centre-line mean). README Step 5.
- verify_cases.py: default tolerance for six written digits set to 1e-5
  (rounding of the written Ceq and T).
- README files: description of the verification conditions corrected
  (five temperatures, four alkali levels, six conditions with fixed p(H2)).

# FACUtility v1.1.0

This release was used to generate all results in the revised paper "FAC-Utility: An Open-Source Framework for Predicting Wall-Thinning Rates in Flow-Accelerated Corrosion".

## Changes

- Newton-Raphson solver for m(H+): the step criterion is now relative to m(H+), and the final step is applied. In the earlier version, the iterations stopped too early, and Ceq was over-predicted by 4-8%.
- Default `newton.absTol` changed from 1e-12 to 1e-24 (code and `FACProperties` of all cases).
- Header of `constant/Kw_vs_T.dat` corrected: the table holds log10(Kw), not ln(Kw).
- Validation case folders renamed to `case1` (2D benchmark), `case2` (3D benchmark), `case3` (3D elbow pipe, previously `DN150`) and `case4`.
- Case 3: `mA` corrected to -1.585e-5 mol/kg.
- Case 4: `mA` corrected to -1.494e-7 mol/kg; `chromiumFactor` set to the calibrated value 0.112.

## Added

- `verification/`: independent check of the chemistry solver (Section 4.7.1 of the paper).
- `postProcessingScripts/`: scripts that make the validation tables.
- `comparisonData/`: measured data and extracted simulation data for Tables 2-11.
- README Section 8: steps to reproduce the tables of the paper.
