#!/bin/bash
# ---------------------------------------------------------------------------
# Phase 5 - Tensor fitting across the cohort.
#
# For every subject: run dtifit on the eddy-corrected data using the ROTATED
# bvecs, then derive AD and RD (dtifit does not write them), then build a
# colour-FA volume for the orientation check.
#
# THE ONE THING THAT MATTERS MOST:
#   dtifit must receive eddy_corrected.eddy_rotated_bvecs, NEVER the original
#   bvec file. Eddy rotated the images to correct motion; the original
#   gradient directions now describe a coordinate frame that no longer
#   exists. Using them produces a normal-LOOKING FA map that is silently
#   wrong. This script refuses to run if the rotated bvecs are absent.
#
# Usage:
#   chmod +x code/phase5_dtifit.sh
#   ./code/phase5_dtifit.sh
#
# Resumable: subjects with an existing dti_FA.nii.gz are skipped.
# Runtime: roughly 1-3 min per subject, so under 2 hours for 78.
# ---------------------------------------------------------------------------

set -u   # not -e: we want to log a failing subject and carry on

PROJ="$HOME/VeronaSemester/DmriMiniCourse"
COHORT="$PROJ/cohort_final.tsv"
PREPROC="$PROJ/derivatives/preproc"
DTIDIR="$PROJ/derivatives/dti"
MOTION_QC="$PROJ/derivatives/motion_qc.tsv"
LOG="$PROJ/logs/phase5_dtifit.log"

EDDY_STEM="eddy_corrected"
MASK_NAME="nodif_brain_mask.nii.gz"     # BET output from Phase 4
RAW="$PROJ/raw"                          # bvals live in the raw BIDS tree

mkdir -p "$DTIDIR" "$(dirname "$LOG")"

# --------------------------- preflight -------------------------------------
echo "=== PREFLIGHT ===" | tee "$LOG"

command -v dtifit >/dev/null 2>&1 || {
  echo "FATAL: dtifit not on PATH. Try:"
  echo "  export PATH=\$PATH:/Users/fredericklancia/fsl/share/fsl/bin"
  exit 1; }

[ -f "$COHORT" ] || { echo "FATAL: no $COHORT"; exit 1; }

if [ ! -f "$MOTION_QC" ]; then
  echo "FATAL: $MOTION_QC not found."
  echo "Run code/phase4f_motion_qc.py first. Fitting tensors for subjects"
  echo "you are about to exclude wastes time and invites you to peek at"
  echo "their FA, which contaminates a pre-specified exclusion rule."
  exit 1
fi

# Subject list = cohort minus motion-QC exclusions.
SUBJECTS=$(awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++){ if($i=="participant_id") id=i;
                              if($i=="exclude_final") ex=i } ; next }
  $ex=="False" || $ex=="0" || $ex=="" { print $id }
' "$MOTION_QC")

N=$(echo "$SUBJECTS" | grep -c . )
echo "Subjects passing motion QC: $N" | tee -a "$LOG"
[ "$N" -gt 0 ] || { echo "FATAL: subject list is empty."; exit 1; }

# Check one subject's files exist before committing to the whole loop.
FIRST=$(echo "$SUBJECTS" | head -1)
resolve_bval () {
  local sub="$1"
  for c in "$RAW/$sub/dwi/${sub}_dwi.bval" \
           "$RAW/$sub/${sub}_dwi.bval" \
           "$PREPROC/$sub/${sub}_dwi.bval" \
           "$PREPROC/$sub/bvals" \
           "$PREPROC/$sub/bvals.txt"; do
    [ -f "$c" ] && { echo "$c"; return 0; }
  done
  for c in "$RAW/$sub"/dwi/*.bval "$RAW/$sub"/*.bval; do
    [ -f "$c" ] && { echo "$c"; return 0; }
  done
  return 1
}

FIRSTBVAL=$(resolve_bval "$FIRST") || {
  echo "FATAL: no .bval file found for $FIRST. Searched under"
  echo "  $RAW/$FIRST/  and  $PREPROC/$FIRST/"
  exit 1; }
echo "bvals resolved to: $FIRSTBVAL"

for f in "$EDDY_STEM.nii.gz" "$EDDY_STEM.eddy_rotated_bvecs" \
         "$MASK_NAME"; do
  if [ ! -f "$PREPROC/$FIRST/$f" ]; then
    echo "FATAL: expected $PREPROC/$FIRST/$f"
    echo "Actual contents of that directory:"
    ls -1 "$PREPROC/$FIRST" 2>/dev/null | sed 's/^/    /'
    echo "Edit the *_NAME variables at the top of this script to match."
    exit 1
  fi
done
echo "Preflight OK on $FIRST" | tee -a "$LOG"

# --------------------------- main loop -------------------------------------
ok=0; skip=0; fail=0
for sub in $SUBJECTS; do
  [[ "$sub" == sub-* ]] || sub="sub-$sub"
  IN="$PREPROC/$sub"
  OUT="$DTIDIR/$sub"
  mkdir -p "$OUT"

  if [ -f "$OUT/dti_RD.nii.gz" ]; then
    skip=$((skip+1)); continue
  fi

  BVEC="$IN/$EDDY_STEM.eddy_rotated_bvecs"
  if [ ! -f "$BVEC" ]; then
    echo "FAIL $sub: no rotated bvecs - refusing to substitute originals" \
      | tee -a "$LOG"
    fail=$((fail+1)); continue
  fi

  # Sanity: rotated bvecs must have the same number of columns as volumes.
  NVOL=$(fslval "$IN/$EDDY_STEM.nii.gz" dim4 | tr -d ' ')
  NBVEC=$(head -1 "$BVEC" | wc -w | tr -d ' ')
  if [ "$NVOL" != "$NBVEC" ]; then
    echo "FAIL $sub: $NVOL volumes but $NBVEC bvec columns" | tee -a "$LOG"
    fail=$((fail+1)); continue
  fi

  BVAL=$(resolve_bval "$sub") || {
    echo "FAIL $sub: no .bval file found" | tee -a "$LOG"
    fail=$((fail+1)); continue; }
  NBVAL=$(tr -s ' \n' ' ' < "$BVAL" | wc -w | tr -d ' ')
  if [ "$NVOL" != "$NBVAL" ]; then
    echo "FAIL $sub: $NVOL volumes but $NBVAL bvals" | tee -a "$LOG"
    fail=$((fail+1)); continue
  fi

  echo "[$(date +%H:%M:%S)] $sub  ($NVOL volumes)" | tee -a "$LOG"

  dtifit \
    -k "$IN/$EDDY_STEM.nii.gz" \
    -o "$OUT/dti" \
    -m "$IN/$MASK_NAME" \
    -r "$BVEC" \
    -b "$BVAL" \
    --save_tensor >> "$LOG" 2>&1

  if [ ! -f "$OUT/dti_FA.nii.gz" ]; then
    echo "FAIL $sub: dtifit produced no FA - see $LOG" | tee -a "$LOG"
    fail=$((fail+1)); continue
  fi

  # AD = L1 (largest eigenvalue, along the fibre).
  # RD = (L2 + L3)/2 (average of the two across-fibre eigenvalues).
  # dtifit writes neither, only L1/L2/L3.
  fslmaths "$OUT/dti_L1" "$OUT/dti_AD"
  fslmaths "$OUT/dti_L2" -add "$OUT/dti_L3" -div 2 "$OUT/dti_RD"

  # Colour-FA volume: |V1| modulated by FA. Channel order is x,y,z, so in
  # standard orientation red=left-right, green=anterior-posterior,
  # blue=superior-inferior.
  fslmaths "$OUT/dti_V1" -abs -mul "$OUT/dti_FA" "$OUT/dti_cFA"

  ok=$((ok+1))
done

# --------------------------- summary ---------------------------------------
echo "" | tee -a "$LOG"
echo "=== PHASE 5 COMPLETE ===" | tee -a "$LOG"
echo "  fitted:  $ok" | tee -a "$LOG"
echo "  skipped: $skip (already had output)" | tee -a "$LOG"
echo "  failed:  $fail" | tee -a "$LOG"

echo ""
echo "Sanity check on the first subject - read these numbers, do not skim:"
S="$DTIDIR/$FIRST"
echo "  FA range (expect ~0 to ~0.9, never >1):"
fslstats "$S/dti_FA" -R | sed 's/^/    /'
echo "  MD mean (expect ~0.0007-0.0010 mm^2/s in brain):"
fslstats "$S/dti_MD" -M | sed 's/^/    /'
echo "  AD mean (should be LARGER than RD mean):"
fslstats "$S/dti_AD" -M | sed 's/^/    /'
echo "  RD mean:"
fslstats "$S/dti_RD" -M | sed 's/^/    /'
echo ""
echo "If FA exceeds 1, or MD is orders of magnitude off, or RD >= AD,"
echo "STOP. Those indicate a units problem, a bad mask, or mis-sorted"
echo "eigenvalues - not something to correct downstream."
echo ""
echo "NEXT: python code/phase5_colorfa_check.py"
