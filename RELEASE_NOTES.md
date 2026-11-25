# FACUtility v1.1.1

No change to the FACUtility code or the validation cases. All results are the same as in v1.1.0.

- Added the owner and boundary files of the FACUtility patch for Cases 1-4
  (comparisonData/simulated/nFail_verification). With these files,
  verify_cases.py checks each face at the temperature of its owner cell (Table 3).
- Added verification/extract_patch_owner.py, which makes these files from the mesh.
- verify_cases.py: default tolerance for six written digits set to 1e-5
  (rounding of the written Ceq and T).
- README files: description of the verification conditions corrected
  (five temperatures, four alkali levels, six conditions with fixed p(H2)).