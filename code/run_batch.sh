#!/bin/bash
# Usage: ./code/run_batch.sh [n_parallel]
set -u
NPAR=${1:-6}
mkdir -p logs

# Iterate over the cohort, not raw/sub-* : raw/ holds ~110 subjects.
subs=$(tail -n +2 cohort_final.tsv | cut -f1)
total=$(echo "$subs" | wc -w | tr -d ' ')
echo "### batch start $(date)  |  ${total} subjects  |  ${NPAR} parallel"

# Only subjects not already complete (resumability guard).
todo=""
for sub in $subs; do
    [ -f "derivatives/preproc/${sub}/eddy_corrected.eddy_rotated_bvecs" ] && continue
    todo="${todo} ${sub}"
done
n_todo=$(echo "$todo" | wc -w | tr -d ' ')
echo "### ${n_todo} to process, $((total - n_todo)) already complete"

echo "$todo" | tr ' ' '\n' | grep -v '^$' | \
    xargs -P "$NPAR" -I{} bash code/process_one_subject.sh {}

echo "### batch end $(date)"

ok=0; bad=""
for sub in $subs; do
    if [ -f "derivatives/preproc/${sub}/eddy_corrected.eddy_rotated_bvecs" ]; then
        ok=$((ok + 1))
    else
        bad="${bad} ${sub}"
    fi
done
echo "complete: ${ok}/${total}"
[ -n "$bad" ] && echo "FAILED:${bad}"
