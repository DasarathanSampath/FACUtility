# Simulated data used in the paper

All values are from FACUtility, with the release named in Section 5.1 of the paper. Cases 1 and 2 are taken at iteration 4000.

Each CSV file is a ParaView export of cell data, one row per wall face, with the face-centre coordinates in `Points:0`, `Points:1` and `Points:2` (see Section 8 of the main README).

| File | Content | Used for |
| --- | --- | --- |
| `case1_meshA.csv` ... `case1_meshF.csv` | Case 1, `downstream` patch, Meshes A-F | Tables 4-6 |
| `case1_meshF_TL_traction_vs_gradient.csv` | Case 1, Mesh F: Uτ from both methods and TL (gradient) | Section 4.6.1 |
| `case2_inflation05.csv`, `case2_inflation08.csv` | Case 2, `downstream` patch, two inflation settings | Tables 4 and 7 |
| `case3_slices/slice_A.csv` ... `slice_E.csv` | Case 3, wall values on the five measuring planes | Tables 8-10 |
| `case4_walls_cr0.112.csv` | Case 4, elbow wall, chromium factor 0.112 | Table 11 |
| `nFail_verification/case1` ... `case4` | `CeqFace`, `T`, `FACProperties`, and `owner`/`boundary` for the patch faces only, of each case | Table 3 |
| `nFail_verification/case1_a0.1`, `case1_a0.01` | the same files for two Case 1 runs with p(H2) fixed at 1e-3 and 1e-6 atm (a = 0.1 and 0.01), plus `pH2_vs_T.dat` of each run | Table 3 |
