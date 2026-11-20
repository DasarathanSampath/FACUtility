# FACUtility [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

FACUtility is an OpenFOAM-based post-processing utility for predicting flow-accelerated corrosion (FAC) wall-thinning rates from a converged CFD solution. The utility combines hydrodynamic quantities derived from the wall region with temperature-dependent chemistry to compute local thinning-related outputs on selected wall patches.

This repository is written for **OpenFOAM V2212** and contains:

- `FACUtility` source code
- required header/source files and build configuration
- four validation / example cases
  - **Case 1: 2D benchmark case** (validationCases/case1)
  - **Case 2: 3D benchmark case** (validationCases/case2)
  - **Case 3: 3D elbow pipe case (plant-representative)** (validationCases/case3)
  - **Case 4: Utanohara & Murase compound elbow** (validationCases/case4)
- the verification test, the comparison data and the scripts that make the tables of the paper (see [Section 8](#8-reproducing-the-tables-of-the-paper))

The default case settings are configured to reproduce the validation conditions reported in the associated publication.

---

## 1. Repository Structure

A typical repository layout is shown below. Adjust the paths if your local folder names differ.

```text
├── FACUtility                     source code (models/, utils/, Make/)
├── validationCases
│   ├── case1                      2D benchmark
│   ├── case2                      3D benchmark
│   ├── case3                      3D elbow pipe
│   └── case4                      Utanohara & Murase compound elbow
├── verification                   independent check of the chemistry solver
├── postProcessingScripts          scripts that make the tables of the paper
├── comparisonData
│   ├── measured                   measured values from the published papers
│   └── simulated                  values extracted from the FACUtility runs
├── CITATION.cff
├── LICENSE
├── README.md
└── RELEASE_NOTES.md
```
---

## 2. Requirements

- **OpenFOAM V2212**
- Linux environment or WSL/Linux-compatible shell with OpenFOAM sourced
- C++ compiler supported by the OpenFOAM installation
- `wmake` available in the environment

Before compiling or running any case, load OpenFOAM V2212:

```bash
source /path/to/OpenFOAM-v2212/etc/bashrc
```

Replace the path above with the correct installation path on your system.

---

## 3. Compiling FACUtility

FACUtility is a standalone custom OpenFOAM utility and is compiled using `wmake`.

### Step 1: Place the source folder

Copy or link the utility source into your user application area, for example:

```bash
mkdir -p $WM_PROJECT_USER_DIR/applications/utilities
cp -r FACUtility $WM_PROJECT_USER_DIR/applications/utilities/
```

### Step 2: Compile

Move into the utility source directory and run `wmake`:

```bash
cd $WM_PROJECT_USER_DIR/applications/utilities/FACUtility
wmake
```

If compilation is successful, the executable will be created in:

```bash
$FOAM_USER_APPBIN/FACUtility
```

### Step 3: Confirm the executable

You can verify that the utility is available using:

```bash
which FACUtility
```

or

```bash
ls $FOAM_USER_APPBIN/FACUtility
```

---

## 4. License

FACUtility is released under the [GNU General Public License v3.0 (GPL-3.0)](LICENSE).

FACUtility directly includes OpenFOAM headers and links against OpenFOAM's core libraries (`fvMesh`, `volScalarField`, `surfaceScalarField`, `fvc::` operators, etc.). Since OpenFOAM itself is licensed under GPL-3.0, FACUtility is a combined/derivative work under the GPL and inherits the same license. See the [`LICENSE`](LICENSE) file for the full text.

## 5. General Workflow

FACUtility is a **decoupled post-processing utility**. 

The workflow is:

1. Set up and run the CFD case using a suitable OpenFOAM solver.
2. Ensure that the flow and temperature fields are converged.
3. Configure `FACProperties` and select the target patch.
4. Run `FACUtility` for the selected time directory.
5. Review the generated post-processing files and ParaView-compatible patch fields.

In summary:

```text
Create mesh -> Run CFD solver -> Run FACUtility -> Post-process outputs
```

---

## 6. Input Requirements

Before running FACUtility, the CFD case must already contain the required solution fields. Depending on the selected physical model and case setup, these typically include:

- velocity field `U`
- pressure field `p` or `p_rgh`
- temperature field `T`
- turbulence quantities required by the parent CFD solver
- mesh and boundary conditions
- `FACProperties` dictionary

The utility reads the converged OpenFOAM case and computes FAC-related quantities only on the specified wall patch.

---

## 7. Running the Four Validation Cases

### 7.1 Case 1: 2D Benchmark

```bash
# move to the case directory
cd validationCases/case1 
# copy the original zero folder as "0"
cp -r 0.orig/ 0 
# create mesh
blockMesh 
# run solver
buoyantSimpleFoam # run solver in single core

# Run in parallel
decomposePar # split mesh, adjust the number of cores in "decomposeParDict"
mpirun -np 4 buoyantSimpleFoam -parallel # run solver in parallel 
mpirun -np 4 buoyantSimpleFoam -parallel > log.solver # to log the solver in a file

mv $(foamListTimes -latestTime)/turbulence*muEff $(foamListTimes -latestTime)/muEff # remove the prefix text from muEff under the latest time directory
FACUtility -latestTime # run FACUtility
gnuplot plots/*.xy # run gnuplots to check the results
```

This case reproduces the 2D validation setup and is intended to demonstrate the local flow acceleration, separation behavior, and resulting FAC-related post-processing.

### 7.2 Case 2: 3D Benchmark

```bash
# move to the case directory
cd validationCases/case2 
# copy the original zero folder as "0"
cp -r 0.orig/ 0 
# create mesh
blockMesh 
# run solver
buoyantSimpleFoam # run solver in single core

# Run in parallel
decomposePar # split mesh, adjust the number of cores in "decomposeParDict"
mpirun -np 4 buoyantSimpleFoam -parallel # run solver in parallel 
mpirun -np 4 buoyantSimpleFoam -parallel > log.solver # to log the solver in a file

mv $(foamListTimes -latestTime)/turbulence*muEff $(foamListTimes -latestTime)/muEff # remove the prefix text from muEff under the latest time directory
FACUtility -latestTime # run FACUtility
gnuplot plots/*.xy # run gnuplots to check the results
```

This case extends the benchmark workflow to a three-dimensional geometry.

### 7.3 Case 3: 3D Elbow Pipe

```bash
# move to the case directory
cd validationCases/case3
# copy the original zero folder as "0"
cp -r 0.orig/ 0 
# create back ground mesh
blockMesh 
# generate *.eMsh for edges
surfaceFeatureExtract
# create mesh
snappyHexMesh -overwrite
# run solver
buoyantSimpleFoam # run solver in single core

# Run in parallel
decomposePar # split mesh, adjust the number of cores in "decomposeParDict"
mpirun -np 4 buoyantSimpleFoam -parallel # run solver in parallel 
mpirun -np 4 buoyantSimpleFoam -parallel > log.solver # to log the solver in a file

mv $(foamListTimes -latestTime)/turbulence*muEff $(foamListTimes -latestTime)/muEff # remove the prefix text from muEff under the latest time directory

FACUtility -latestTime # run FACUtility

Note: The plots are visualized and generated using ParaView
```

This case reproduces the elbow-pipe validation configuration and is useful for studying circumferential variation in thinning around the bend.

### 7.4 Case 4: Utanohara and Murase Compound Elbow

```bash
# move to the case directory
cd validationCases/case4
# copy the original zero folder as "0"
cp -r 0.orig/ 0 
# create back ground mesh
blockMesh 
# generate *.eMsh for edges
surfaceFeatureExtract
# create mesh
snappyHexMesh -overwrite
# run solver
buoyantSimpleFoam # run solver in single core

# Run in parallel
decomposePar # split mesh, adjust the number of cores in "decomposeParDict"
mpirun -np 4 buoyantSimpleFoam -parallel # run solver in parallel 
mpirun -np 4 buoyantSimpleFoam -parallel > log.solver # to log the solver in a file

mv $(foamListTimes -latestTime)/turbulence*muEff $(foamListTimes -latestTime)/muEff # remove the prefix text from muEff under the latest time directory

FACUtility -latestTime # run FACUtility

Note: The plots are visualized and generated using ParaView
```


---

## 8. Reproducing the Tables of the Paper

These steps re-create the verification test and the validation tables of the paper from the data in `comparisonData/`. They need only Python. OpenFOAM is not needed.

### Step 1: Install Python packages

```bash
pip install numpy pandas mpmath
```

Python 3.8 or newer is required. Run all commands below from the repository root folder.

```bash
mkdir -p results
```

### Step 2: Chemistry verification (Tables 2 and 3)

```bash
# independent solution of Eqs. (6)-(9) at 21 conditions (Table 2)
python3 verification/verify_chemistry.py --out results/reference_values.csv

# self-test: the first must PASS 4 of 4, the second must FAIL 4 of 4
python3 verification/verify_chemistry.py --compare verification/selftest_correct.csv
python3 verification/verify_chemistry.py --compare verification/selftest_k2div.csv

# check of CeqFace from the four validation runs, and failed solves (Table 3)
python3 verification/verify_cases.py --root comparisonData/simulated/nFail_verification
```

### Step 3: Case 1, mesh sensitivity (Tables 4, 5 and 6)

```bash
python3 postProcessingScripts/mesh_window_stats.py \
    comparisonData/simulated/case1_meshA.csv comparisonData/simulated/case1_meshB.csv \
    comparisonData/simulated/case1_meshC.csv comparisonData/simulated/case1_meshD.csv \
    comparisonData/simulated/case1_meshE.csv comparisonData/simulated/case1_meshF.csv \
    --origin 0.253 --out results/case1_mesh_window_stats.csv
```

### Step 4: Case 1, traction and gradient methods (Section 4.6.1)

```bash
python3 postProcessingScripts/case1_traction_vs_gradient.py \
    comparisonData/simulated/case1_meshF_TL_traction_vs_gradient.csv \
    --out results/case1_meshF_TL_paired_window.csv
```

### Step 5: Case 2 (Tables 4 and 7)

```bash
python3 postProcessingScripts/mesh_window_stats.py \
    comparisonData/simulated/case2_inflation05.csv comparisonData/simulated/case2_inflation08.csv \
    --span-coord Points:2 --origin 0.253 --out results/case2_mesh_window_stats.csv
```

### Step 6: Case 3 (Tables 8, 9 and 10)

```bash
python3 postProcessingScripts/process_elbow_slices_4.py comparisonData/simulated/case3_slices
python3 postProcessingScripts/case3_compare.py
```

### Step 7: Case 4 (Table 11)

```bash
python3 postProcessingScripts/case4_sensors.py \
    comparisonData/simulated/case4_walls_cr0.112.csv --factor 0.112 \
    --out results/case4_sensor_values.csv
```

### Expected results

| Step | What to check | Expected value |
| --- | --- | --- |
| 2 | `verify_cases.py` summary | case1 PASS, case2 PASS, case3 PASS, case4 PASS; 0 fallback faces |
| 3 | Mesh F, mean TL over 5-20 mm | 0.7910 mm/10,000 h |
| 4 | Window-mean TL, gradient -> traction | 0.7910 -> 0.7911 (+0.005 %) |
| 5 | Inflation 0.5, centre-line mean TL | 0.643 mm/10,000 h |
| 6 | RMSE / r² of the normalized comparison | 0.177 / 0.454 |
| 7 | Mean sim/exp, eight straight-pipe sensors | 0.99 |

---

### Notes on the steps

**Step 2.** `verify_chemistry.py` solves the chemistry model directly from Eqs. (6)-(9) of the paper. It does not use any FACUtility code. The hydrogen-ion balance is solved with 50-digit precision by two methods, and the two must agree. Two of the 21 conditions use p(H2) = 10⁻³ and 10⁻⁶ atm (a = 0.1 and 0.01). At these conditions, multiplying and dividing the Fe(OH)2 term by `a` give very different results, so the test can tell them apart. `selftest_k2div.csv` holds values made with the wrong form (`k2/a`). It must fail, which shows that the test can detect this error. `verify_cases.py` checks every face of `CeqFace` against the independent solution at the cell temperature. It also counts faces with the fallback value of a failed solve. More details are in `verification/README.md`.

**Steps 3 and 5.** The comparison window is 5-20 mm downstream of the upstream edge of the specimen, which is at x = 0.253 m in Cases 1 and 2 (`--origin`). For each mesh, the script gives the mean, minimum and maximum TL in the window, the y+ distribution, and the share of faces with 30 ≤ y+ ≤ 50. For Case 2, `--span-coord` gives the mean over a 2 mm wide strip along the centre line, which is where the experiment measured.

**Step 4.** The thinning rate from the traction-based shear velocity is found face by face from Eq. (15): TL_t = TL_g (Uτ,t / Uτ,g)². All means are arithmetic means over the faces, which have the same size.

**Step 6.** `process_elbow_slices_4.py` reads the five slice files and finds TL at the 40 measuring points (5 planes × 8 angles). `case3_compare.py` normalizes these values by their maximum, compares them with the measured values, and writes the tables to `processed_slices/`.

**Step 7.** `--factor` is the chromium factor used in the run that produced the CSV file. Each sensor value is the mean TL over the wall faces in a 20 mm × 3 mm area at the sensor position. The geometry (R = 76.2 mm, straight section 54 mm, D = 49.5 mm) is the same as the test section. The script also lists the faces with y+ > 70 and checks that none of them is inside a sensor area.

**ParaView export.** All CSV files in `comparisonData/simulated/` were exported from ParaView as cell data, one row per wall face: tick only the wall patch in the reader, apply **Filters > Cell Centers** (with *Vertex Cells* on), then **File > Save Data > CSV**. The columns `Points:0`, `Points:1` and `Points:2` hold the face-centre coordinates.
