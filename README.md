# ColliderBit script for Master’s thesis

This repository contains scripts and analysis files used with **GAMBIT / ColliderBit** for a Master’s thesis project studying collider constraints on simplified models and MSSM parameter scans.

---

## Structure

The `Data/` directory is organised into two main parts:

- `m_vs_BR/`  
  Contains results for the simplified model analysis (mass vs branching ratio).

- `mu_vs_M2/`  
  Contains results for the MSSM parameter scan (μ vs M₂ plane).

The corresponding SLHA spectrum files and YAML configuration files used in the runs are stored in their respective folders.

---

## Prerequisites

This project requires **GAMBIT (ColliderBit module)**.

Install GAMBIT following the official documentation:

- https://gambitbsm.org/documentation/installation/introduction/
- https://github.com/GambitBSM/gambit_2.5

Make sure ColliderBit is enabled during installation and compilation.

---

## Setup

### 1. Add analysis code

Place the custom ColliderBit analysis C++ file into:



Register the analysis in the ColliderBit analysis container in the same directory, following the GAMBIT analysis framework structure.

---

### 2. Build GAMBIT

Build GAMBIT with ColliderBit enabled using the official instructions:

https://gambitbsm.org/documentation/installation/building/

---

### 3. Input files

Ensure the following files are placed in the run directory:

- YAML file (scan configuration)
- SLHA files (model spectra)

These define the parameter scan and physics setup used in ColliderBit.

---

### 4. Run

Run GAMBIT with:

```bash
./gambit input.yaml
