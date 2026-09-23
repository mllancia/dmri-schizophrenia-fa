#!/bin/bash
# ---------------------------------------------------------------------------
# Phase 8 - TBSS voxelwise analysis.
#
# THE ONE THING THAT MATTERS MOST:
#   randomise matches design-matrix rows to subjects purely by POSITION.
#   TBSS merges subjects in alphabetical filename order. If the design
#   matrix is built in any other order, diagnosis labels get attached to the
#   wrong brains and the result is meaningless - with no error message.
#   This script builds both from the same sorted list and verifies they agree.
#
# Model: FA ~ CONTROL + SCZ + age + sex + motion   (covariates demeaned)
# Contrasts: 1 = CONTROL > SCZ  (the hypothesised direction)
#            2 = SCZ > CONTROL
#
# Usage:
#   chmod +x code/phase8_tbss.sh
#   caffeinate -i nohup ./code/phase8_tbss.sh > logs/phase8.log 2>&1 &
#   tail -f logs/phase8.log
#
# Runtime: tbss_2_reg is the slow step on a cluster but was fast for fnirt
# here; randomise with 5000 permutations + TFCE is likely the longest step.
# ---------------------------------------------------------------------------

set -u
PROJ="$HOME/VeronaSemester/DmriMiniCourse"
cd "$PROJ" || exit 1
TB="$PROJ/derivatives/tbss"
NPERM="${1:-5000}"

command -v tbss_1_preproc >/dev/null || { echo "FATAL: TBSS not on PATH"; exit 1; }
command -v randomise      >/dev/null || { echo "FATAL: randomise not on PATH"; exit 1; }

if [ -d "$TB/stats" ] && [ -f "$TB/stats/all_FA_skeletonised.nii.gz" ]; then
  echo "TBSS stages 1-4 already complete; skipping to randomise."
else
  rm -rf "$TB"; mkdir -p "$TB"

  # ---------------- inputs: analysis cohort ONLY ----------------
  source "$PROJ/.venv/bin/activate"
  python - <<'EOF'
import pandas as pd, shutil, os, numpy as np
P = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
q = pd.read_csv(f"{P}/derivatives/motion_qc.tsv", sep="\t")
q["exclude_final"] = q.exclude_final.astype(bool)
k = q.loc[~q.exclude_final].copy()
c = pd.read_csv(f"{P}/cohort_final.tsv", sep="\t")
k = k.merge(c[["participant_id","age","gender"]], on="participant_id")
k = k.sort_values("participant_id").reset_index(drop=True)
assert len(k) == 72, f"expected 72, got {len(k)}"

for s in k.participant_id:
    src = f"{P}/derivatives/dti/{s}/dti_FA.nii.gz"
    assert os.path.isfile(src), src
    shutil.copy(src, f"{P}/derivatives/tbss/{s}_FA.nii.gz")

con = (k.diagnosis == "CONTROL").astype(int)
scz = (k.diagnosis == "SCHZ").astype(int)
age = k.age - k.age.mean()
sex = (k.gender.astype(str).str.upper().str[0] == "M").astype(float)
sex = sex - sex.mean()
mot = k.mean_rel_rms - k.mean_rel_rms.mean()
X = np.column_stack([con, scz, age, sex, mot])

with open(f"{P}/derivatives/tbss/design.mat", "w") as f:
    f.write(f"/NumWaves 5\n/NumPoints {len(k)}\n/Matrix\n")
    for r in X: f.write(" ".join(f"{v:.6f}" for v in r) + "\n")
with open(f"{P}/derivatives/tbss/design.con", "w") as f:
    f.write("/ContrastName1 CONTROL>SCZ\n/ContrastName2 SCZ>CONTROL\n"
            "/NumWaves 5\n/NumContrasts 2\n/Matrix\n"
            "1 -1 0 0 0\n-1 1 0 0 0\n")
k[["participant_id","diagnosis","age","gender","mean_rel_rms"]].to_csv(
    f"{P}/derivatives/tbss/design_subjects.tsv", sep="\t", index=False)
print(f"copied {len(k)} FA maps; design {X.shape}; "
      f"{int(con.sum())} CONTROL / {int(scz.sum())} SCZ")
EOF
  [ $? -eq 0 ] || { echo "FATAL: input preparation failed"; exit 1; }

  cd "$TB"
  echo "[$(date +%H:%M)] tbss_1_preproc";  tbss_1_preproc *.nii.gz
  echo "[$(date +%H:%M)] tbss_2_reg -T";   tbss_2_reg -T
  echo "[$(date +%H:%M)] tbss_3_postreg -S"; tbss_3_postreg -S
  echo "[$(date +%H:%M)] tbss_4_prestats 0.2"; tbss_4_prestats 0.2
  cd "$PROJ"
fi

# ---------------- order verification: do not skip ----------------
echo ""
echo "=== ORDER CHECK ==="
ls -1 "$TB/FA/"*_FA.nii.gz | sed -E 's|.*/||; s|(_FA)+\.nii\.gz$||' > /tmp/tbss_order.txt
tail -n +2 "$TB/design_subjects.tsv" | cut -f1 > /tmp/design_order.txt
if diff -q /tmp/tbss_order.txt /tmp/design_order.txt >/dev/null; then
  echo "OK: TBSS merge order matches design matrix row order ($(wc -l < /tmp/tbss_order.txt | tr -d ' ') subjects)"
else
  echo "FATAL: ORDER MISMATCH between TBSS inputs and design matrix."
  diff /tmp/tbss_order.txt /tmp/design_order.txt | head
  exit 1
fi
NVOL=$(fslval "$TB/stats/all_FA_skeletonised" dim4 | tr -d ' ')
[ "$NVOL" = "72" ] || { echo "FATAL: all_FA_skeletonised has $NVOL volumes"; exit 1; }

# ---------------- randomise ----------------
echo ""
echo "[$(date +%H:%M)] randomise: $NPERM permutations, TFCE"
cd "$TB/stats"
randomise -i all_FA_skeletonised -o tbss_FA \
          -m mean_FA_skeleton_mask \
          -d ../design.mat -t ../design.con \
          -n "$NPERM" --T2 -D

echo ""
echo "=== PHASE 8 COMPLETE [$(date +%H:%M)] ==="
for t in 1 2; do
  f="tbss_FA_tfce_corrp_tstat${t}.nii.gz"
  [ -f "$f" ] || continue
  MAX=$(fslstats "$f" -R | awk '{print $2}')
  NSIG=$(fslstats "$f" -l 0.95 -V | awk '{print $1}')
  NSK=$(fslstats mean_FA_skeleton_mask -V | awk '{print $1}')
  lab=$([ $t = 1 ] && echo "CONTROL>SCZ" || echo "SCZ>CONTROL")
  echo "  contrast $t ($lab): max 1-p = $MAX, significant voxels $NSIG of $NSK"
done
echo ""
echo "Threshold is 1-p > 0.95 (randomise outputs 1-p, not p)."
echo "NEXT: python code/phase8b_tbss_report.py"
