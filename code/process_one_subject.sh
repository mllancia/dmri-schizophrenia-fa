#!/bin/bash
# Preprocess one subject: MP-PCA -> Gibbs -> mask -> eddy
# Usage: ./code/process_one_subject.sh sub-10193
# Run from the project root. Assumes the venv is already active.

set -u
sub=$1
root=$(pwd)
out="derivatives/preproc/${sub}"
log="logs/${sub}.log"
mkdir -p "$out" logs qc

echo "=== ${sub} started $(date '+%H:%M:%S')" | tee -a "$log"

# --- 4a/4b: denoise + de-Gibbs (script has its own skip guard) ---
python code/denoise_one.py "$sub" >> "$log" 2>&1 || {
    echo "FAIL ${sub}: denoise" | tee -a "$log"; exit 1; }

cd "$out" || exit 1

# --- 4c: brain mask from the b=0 volume (index 0, verified) ---
if [ ! -f nodif_brain_mask.nii.gz ]; then
    fslroi "${sub}_degibbs.nii.gz" nodif 0 1 >> "${root}/${log}" 2>&1
    bet nodif nodif_brain -m -f 0.25   >> "${root}/${log}" 2>&1
fi
[ -f nodif_brain_mask.nii.gz ] || {
    echo "FAIL ${sub}: mask" | tee -a "${root}/${log}"; exit 1; }

# --- 4d: eddy (no --topup: single PE direction, no fieldmap) ---
if [ ! -f eddy_corrected.nii.gz ]; then
    eddy_cpu \
      --imain="${sub}_degibbs.nii.gz" \
      --mask=nodif_brain_mask.nii.gz \
      --acqp="${root}/acqp.txt" \
      --index="${root}/index.txt" \
      --bvecs="${root}/raw/${sub}/dwi/${sub}_dwi.bvec" \
      --bvals="${root}/raw/${sub}/dwi/${sub}_dwi.bval" \
      --repol \
      --out=eddy_corrected >> "${root}/${log}" 2>&1
fi

cd "$root" || exit 1

if [ -f "${out}/eddy_corrected.eddy_rotated_bvecs" ]; then
    echo "=== ${sub} OK $(date '+%H:%M:%S')" | tee -a "$log"
else
    echo "FAIL ${sub}: eddy" | tee -a "$log"; exit 1
fi
