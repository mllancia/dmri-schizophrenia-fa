#!/usr/bin/env python3
"""
Phase 5b - Cross-validate FSL's tensor fit against the from-scratch
implementation (validate_fit.py, salvaged from the Brusini course project).

WHY
---
FSL is a black box here. You independently implemented the same least-squares
tensor fit. If the two agree at the same voxel, that is a positive control
demonstrating the pipeline is doing what you think it is, and it is worth a
paragraph in the report. If they disagree, one of them is wrong and you need
to know which before interpreting anything.

The fits will not agree to machine precision - FSL may weight differently or
handle the b=0 rows differently - but FA should agree to roughly two decimal
places in well-conditioned white matter voxels.

Usage:
    source .venv/bin/activate
    python code/phase5b_crossvalidate.py                # picks voxels, prints
    python code/phase5b_crossvalidate.py sub-10321 5    # subject, n voxels
"""

import os
import sys
import subprocess
import numpy as np
import nibabel as nib

PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
os.chdir(PROJ)

SUB = sys.argv[1] if len(sys.argv) > 1 else "sub-10321"
NVOX = int(sys.argv[2]) if len(sys.argv) > 2 else 5

DTI = f"derivatives/dti/{SUB}"
PRE = f"derivatives/preproc/{SUB}"
DWI = f"{PRE}/eddy_corrected.nii.gz"
BVEC = f"{PRE}/eddy_corrected.eddy_rotated_bvecs"
BVAL = f"raw/{SUB}/dwi/{SUB}_dwi.bval"

for f in (DWI, BVEC, BVAL, f"{DTI}/dti_FA.nii.gz"):
    if not os.path.isfile(f):
        sys.exit(f"FATAL: missing {f}")

fa = nib.load(f"{DTI}/dti_FA.nii.gz").get_fdata()
md = nib.load(f"{DTI}/dti_MD.nii.gz").get_fdata()
l1 = nib.load(f"{DTI}/dti_L1.nii.gz").get_fdata()
l2 = nib.load(f"{DTI}/dti_L2.nii.gz").get_fdata()
l3 = nib.load(f"{DTI}/dti_L3.nii.gz").get_fdata()

# Pick well-conditioned white matter voxels: high FA, plausible MD, all
# eigenvalues positive. Deterministic seed so the report is reproducible.
good = (fa > 0.55) & (fa < 0.95) & (md > 0.0005) & (md < 0.0012) \
       & (l1 > 0) & (l2 > 0) & (l3 > 0)
idx = np.argwhere(good)
if len(idx) < NVOX:
    sys.exit(f"FATAL: only {len(idx)} candidate voxels found")
rng = np.random.default_rng(42)
picks = idx[rng.choice(len(idx), NVOX, replace=False)]

print(f"Cross-validation: {SUB}, {NVOX} white matter voxels\n")
print(f"{'voxel (i,j,k)':<18}{'FSL FA':>9}{'FSL MD':>11}"
      f"{'FSL L1':>11}{'FSL L3':>11}")
print("-" * 60)
for i, j, k in picks:
    print(f"{f'({i},{j},{k})':<18}{fa[i,j,k]:>9.4f}{md[i,j,k]:>11.6f}"
          f"{l1[i,j,k]:>11.6f}{l3[i,j,k]:>11.6f}")

print("\nNow run your own implementation at each voxel:\n")
for i, j, k in picks:
    print(f"  python validate_fit.py {DWI} {BVAL} {BVEC} {i} {j} {k}")

print("\nCompare FA and MD. Agreement to ~2 decimal places is expected and")
print("sufficient. A systematic discrepancy (e.g. all your FA values higher)")
print("points at a convention difference - check b-value units, the sign of")
print("the log, and whether both fits include the b=0 rows.")
print("\nIf validate_fit.py takes different argument names or order, adjust")
print("the printed commands to match.")

# Attempt to run it automatically if the script is present.
if os.path.isfile("validate_fit.py"):
    print("\n" + "=" * 60)
    print("Attempting automatic comparison")
    print("=" * 60)
    for i, j, k in picks:
        try:
            r = subprocess.run(
                [sys.executable, "validate_fit.py",
                 DWI, BVAL, BVEC, str(i), str(j), str(k)],
                capture_output=True, text=True, timeout=120)
            out = (r.stdout or r.stderr).strip().replace("\n", " | ")
            print(f"({i},{j},{k}) FSL FA={fa[i,j,k]:.4f} -> {out[:160]}")
        except Exception as e:
            print(f"({i},{j},{k}) could not run: {e}")
            break
else:
    print("\nvalidate_fit.py not found - locate it and rerun, or use "
          "the commands above manually.")
