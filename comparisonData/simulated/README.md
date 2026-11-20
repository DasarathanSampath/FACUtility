# Simulated data used in the paper

All values are from FACUtility at iteration 4000 (Cases 1 and 2), with the release named in Section 5.1 of the paper.

| File | Made with | Paper |
| --- | --- | --- |
| `case1_meshA.csv` ... `case1_meshF.csv` | ParaView export of the `downstream` patch | input to Tables 4-6 |
| `case1_mesh_window_stats.csv` | `mesh_window_stats.py` | Tables 4-6 |
| `case2_inflation05.csv`, `case2_inflation08.csv` | ParaView export of the `downstream` patch | input to Tables 4 and 7 |
| `case2_mesh_window_stats.csv` | `mesh_window_stats.py --span-coord Points:2` | Tables 4 and 7 |
| `case1_meshF_TL_traction_vs_gradient.csv` | paired TL, face by face, in the 5-20 mm window | Section 4.6.1 |
| `case3_slices/` | ParaView export of the five slices | input to Tables 8-9 |
| `case3_comparison_unrounded.csv` | `case3_compare.py` | Tables 8-9 (unrounded) |
| `case3_tables_8_9.csv` | `case3_compare.py` | Tables 8-9 |
| `case3_statistics.csv` | `case3_compare.py` | Table 10 |
| `case4_walls_cr0.112.csv` | ParaView export of the elbow wall | input to Table 11 |
| `case4_sensor_values.csv` | `case4_sensors.py --factor 0.112` | Table 11 |
| `nFail_verification/case1-4/` | `CeqFace`, `T`, `FACProperties` of each case | Table 3 |
| `nFail_verification/verify_cases_output.txt` | `verify_cases.py` | Table 3 |
