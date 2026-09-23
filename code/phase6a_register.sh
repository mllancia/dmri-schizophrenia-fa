#!/bin/bash
# ---------------------------------------------------------------------------
# Phase 6a - Register every subject's diffusion maps into FMRIB58_FA standard
# space, so the JHU atlas can be applied to all of them in a common frame.
#
# ROUTE B: FA -> FMRIB58_FA directly. The subject's T1 is never touched, which
# is what keeps the group-imbalanced ghosting artifact out of the pipeline.
#
# Two stages per subject:
#   flirt  - 12-DOF affine, gets global position/scale/shear approximately right
#   fnirt  - nonlinear warp using FSL's FA_2_FMRIB58_1mm config, which is tuned
#            for exactly this registration
# Then applywarp carries FA, MD, AD and RD through the same warp.
#
# Usage:
#   chmod +x code/phase6a_register.sh
#   caffeinate -i nohup ./code/phase6a_register.sh 5 > logs/phase6a.log 2>&1 &
#   tail -f logs/phase6a.log
#
# Argument: number of parallel jobs (default 5).
# Runtime: ~5-15 min per subject for fnirt, so ~2-4 hours at 5-way for 72.
# Resumable: subjects with an existing FA_std.nii.gz are skipped.
# ---------------------------------------------------------------------------

set -u

PROJ="$HOME/VeronaSemester/DmriMiniCourse"
cd "$PROJ" || exit 1

FSLDIR="${FSLDIR:-$HOME/fsl}"
TARGET="$FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz"
CONFIG="FA_2_FMRIB58_1mm"
DTIDIR="$PROJ/derivatives/dti"
STDDIR="$PROJ/derivatives/std"
LOG="$PROJ/logs/phase6a_register.log"

# ----- worker: called by xargs, one subject per invocation ------------------
if [ "${1:-}" = "--one" ]; then
  sub="$2"
  IN="$DTIDIR/$sub"
  OUT="$STDDIR/$sub"
  mkdir -p "$OUT"

  [ -f "$OUT/FA_std.nii.gz" ] && { echo "skip $sub"; exit 0; }

  echo "[$(date +%H:%M:%S)] start $sub"

  flirt -in "$IN/dti_FA" -ref "$TARGET" \
        -omat "$OUT/fa2std_affine.mat" -dof 12 \
        -out "$OUT/FA_affine" >> "$LOG" 2>&1 || {
    echo "FAIL $sub: flirt"; exit 1; }

  fnirt --in="$IN/dti_FA" --ref="$TARGET" \
        --aff="$OUT/fa2std_affine.mat" \
        --config="$CONFIG" \
        --cout="$OUT/fa2std_warp" \
        --iout="$OUT/FA_std" >> "$LOG" 2>&1 || {
    echo "FAIL $sub: fnirt"; exit 1; }

  for map in MD AD RD; do
    applywarp -i "$IN/dti_$map" -r "$TARGET" \
              -w "$OUT/fa2std_warp" -o "$OUT/${map}_std" >> "$LOG" 2>&1
  done

  if [ -f "$OUT/FA_std.nii.gz" ]; then
    echo "[$(date +%H:%M:%S)] done  $sub"
  else
    echo "FAIL $sub: no FA_std written"
  fi
  exit 0
fi

# ----- driver --------------------------------------------------------------
NPROC="${1:-5}"
mkdir -p "$STDDIR" "$(dirname "$LOG")"

echo "=== PHASE 6a PREFLIGHT ==="
command -v fnirt >/dev/null 2>&1 || { echo "FATAL: fnirt not on PATH"; exit 1; }
[ -f "$TARGET" ] || { echo "FATAL: target not found: $TARGET"; exit 1; }
[ -f "$FSLDIR/etc/flirtsch/$CONFIG.cnf" ] || {
  echo "FATAL: fnirt config not found: $FSLDIR/etc/flirtsch/$CONFIG.cnf"
  echo "Available configs:"; ls -1 "$FSLDIR/etc/flirtsch/"*.cnf | sed 's/^/  /'
  exit 1; }
[ -f "derivatives/motion_qc.tsv" ] || {
  echo "FATAL: derivatives/motion_qc.tsv missing"; exit 1; }

# Subject list from motion_qc.tsv only - NEVER from a directory glob.
SUBJECTS=$(awk -F'\t' '
  NR==1 { for(i=1;i<=NF;i++){ if($i=="participant_id") id=i;
                              if($i=="exclude_final") ex=i } ; next }
  tolower($ex)=="false" || $ex=="0" { print $id }
' derivatives/motion_qc.tsv)

N=$(echo "$SUBJECTS" | grep -c .)
echo "Analysis cohort: $N subjects"
echo "Target:  $TARGET"
echo "Config:  $CONFIG"
echo "Jobs:    $NPROC parallel"
[ "$N" -eq 72 ] || echo "!! expected 72 - check motion_qc.tsv exclusions"
echo ""

# xargs rather than 'wait -n', which macOS bash 3.2 does not support.
echo "$SUBJECTS" | xargs -P "$NPROC" -I{} "$0" --one {}

echo ""
echo "=== PHASE 6a COMPLETE ==="
DONE=$(ls -1d "$STDDIR"/sub-*/FA_std.nii.gz 2>/dev/null | wc -l | tr -d ' ')
echo "FA_std written for $DONE of $N subjects"
echo ""
echo "NOW CHECK THE REGISTRATION - do not skip this:"
echo "  fslmerge -t /tmp/all_FA_std \$(ls -1d $STDDIR/sub-*/FA_std.nii.gz)"
echo "  fslmaths /tmp/all_FA_std -Tmean /tmp/mean_FA_std"
echo "  fsleyes $TARGET /tmp/mean_FA_std -cm hot &"
echo ""
echo "A sharp mean FA map with recognisable tract structure means registration"
echo "worked. A blurred smear means it did not, and every downstream ROI value"
echo "would be averaged across mismatched anatomy."
echo ""
echo "NEXT: python code/phase6b_roi_extract.py"
