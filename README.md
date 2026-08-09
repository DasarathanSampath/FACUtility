# FACUtility [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

FACUtility is an OpenFOAM-based post-processing utility for predicting flow-accelerated corrosion (FAC) wall-thinning rates from a converged CFD solution. The utility combines hydrodynamic quantities derived from the wall region with temperature-dependent chemistry to compute local thinning-related outputs on selected wall patches.

This repository is written for **OpenFOAM V2212** and contains:

- `FACUtility` source code
- required header/source files and build configuration
- four validation / example cases
  - **2D benchmark case** (validationCases/2DCase)
  - **3D benchmark case** (validationCases/3DCase)
  - **3D elbow pipe case (plant-representative)** (validationCases/DN150)
  - **Utanohara & Murase compound elbow** (validationCases/case4)

The default case settings are configured to reproduce the validation conditions reported in the associated publication.

---

## 1. Repository Structure

A typical repository layout is shown below. Adjust the paths if your local folder names differ.

```text
├── FACUtility
│   ├── FACUtility.C
│   ├── Make
│   │   ├── files
│   │   └── options
│   ├── createFields.H
│   ├── models
│   │   ├── FACChemistryConstants.C
│   │   ├── FACChemistryConstants.H
│   │   ├── FACChemistryModel.C
│   │   ├── FACChemistryModel.H
│   │   ├── FACErosionModel.C
│   │   ├── FACErosionModel.H
│   │   ├── FACNewtonRaphson.H
│   │   ├── FACThinningModel.C
│   │   └── FACThinningModel.H
│   ├── readFACProperties.H
│   ├── utils
│   │   ├── FACCSV.H
│   │   ├── FACChemistryTable.C
│   │   ├── FACChemistryTable.H
│   │   ├── FACDatWriter.H
│   │   ├── FACInterpolation.C
│   │   ├── FACInterpolation.H
│   │   ├── FACMapping.H
│   │   ├── FACProperties.H
│   │   ├── FACSelection.C
│   │   └── FACSelection.H
│   └── writeOutputs.H
├── LICENSE
├── README.md
└── validationCases
    ├── 2DCase
    │   ├── 0.orig
    │   ├── constant
    │   ├── plots
    │   └── system
    ├── 3DCase
    │   ├── 0.orig
    │   ├── constant
    │   ├── plots
    │   └── system
    ├── case4
    │   ├── 0.orig
    │   ├── constant
    │   ├── plots
    │   └── system
    └── DN150
        ├── 0.orig
        ├── constant
        └── system
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

### 7.1 2D Benchmark Case

```bash
# move to the case director
cd validationCases/2DCase 
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

### 7.2 3D Benchmark Case

```bash
# move to the case director
cd validationCases/3DCase 
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

### 7.3 DN150

```bash
# move to the case director
cd validationCases/DN150
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

### 7.4 Case 4 (Utanohara and Murase)

```bash
# move to the case director
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

