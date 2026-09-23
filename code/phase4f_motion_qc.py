#!/usr/bin/env python3
"""
Phase 4f - Motion QC across the 78-subject cohort.

Reads eddy's motion outputs for every subject in cohort_final.tsv, computes
three motion metrics per subject, applies PRE-SPECIFIED exclusion thresholds,
and tests whether motion differs by diagnosis.

WHY THIS EXISTS
---------------
People with schizophrenia move more in scanners on average. Motion lowers
measured FA. So motion alone can manufacture a perfect "reduced FA in SCZ"
result with no biology in it. This script (a) removes the worst subjects and
(b) produces the covariate you will put in the Phase 7 GLM.

RUN BEFORE dtifit. Do not look at any FA value until this has been run and
the thresholds below have been frozen.

Usage:
    source ~/VeronaSemester/DmriMiniCourse/.venv/bin/activate
    python code/phase4f_motion_qc.py
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# CONFIGURATION - verify these against your actual tree before running
# ---------------------------------------------------------------------------
PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
COHORT = os.path.join(PROJ, "cohort_final.tsv")
PREPROC = os.path.join(PROJ, "derivatives", "preproc")   # <- subject dirs live here
EDDY_STEM = "eddy_corrected"                             # <- eddy's -out prefix
OUT_TSV = os.path.join(PROJ, "derivatives", "motion_qc.tsv")

# ---------------------------------------------------------------------------
# PRE-SPECIFIED EXCLUSION THRESHOLDS
# ---------------------------------------------------------------------------
# CHECK DTI_todo_annotated.md FIRST. If you already wrote thresholds there,
# use those numbers, not these. Changing thresholds after seeing FA is
# p-hacking. These defaults follow common practice in the DTI literature.
THRESH_MEAN_REL_RMS = 2.0    # mm, mean volume-to-volume displacement
THRESH_MAX_REL_RMS = 3.0     # mm, worst single volume-to-volume jump
THRESH_PCT_OUTLIER = 10.0    # %, slices eddy flagged and replaced

# Matched-pair design: if a subject is dropped, drop its partner too?
DROP_PARTNER = True


def read_movement_rms(path):
    """eddy's .eddy_movement_rms: one row per volume, two columns.
    col 0 = RMS displacement relative to the FIRST volume (absolute)
    col 1 = RMS displacement relative to the PREVIOUS volume (relative)
    Relative is the one that tracks jerky motion, so it is our primary metric.
    """
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return arr[:, 0], arr[:, 1]


def read_outlier_map(path):
    """eddy's .eddy_outlier_map: one header line, then a volumes x slices
    matrix of 0/1. A 1 means eddy judged that slice an outlier and replaced
    it (because --repol was on)."""
    with open(path) as fh:
        lines = fh.readlines()
    body = [l for l in lines[1:] if l.strip()]
    return np.array([[int(v) for v in l.split()] for l in body])


def main():
    if not os.path.isfile(COHORT):
        sys.exit(f"FATAL: cohort file not found: {COHORT}")

    cohort = pd.read_csv(COHORT, sep="\t")

    # Be tolerant about column naming across earlier scripts.
    id_col = next((c for c in ("participant_id", "subject", "sub", "id")
                   if c in cohort.columns), None)
    dx_col = next((c for c in ("diagnosis", "group", "dx")
                   if c in cohort.columns), None)
    pair_col = next((c for c in ("pair_id", "pair", "match_id")
                     if c in cohort.columns), None)

    if id_col is None or dx_col is None:
        sys.exit(f"FATAL: could not find id/diagnosis columns. "
                 f"Found: {list(cohort.columns)}")

    print(f"Cohort: {len(cohort)} subjects  "
          f"(id={id_col}, dx={dx_col}, pair={pair_col})")

    rows, missing = [], []

    for _, r in cohort.iterrows():
        sub = str(r[id_col])
        if not sub.startswith("sub-"):
            sub = "sub-" + sub
        d = os.path.join(PREPROC, sub)

        f_rms = os.path.join(d, EDDY_STEM + ".eddy_movement_rms")
        f_out = os.path.join(d, EDDY_STEM + ".eddy_outlier_map")

        if not os.path.isfile(f_rms):
            missing.append((sub, f_rms))
            continue

        try:
            abs_rms, rel_rms = read_movement_rms(f_rms)
        except Exception as e:
            missing.append((sub, f"unreadable {f_rms}: {e}"))
            continue

        pct_out = np.nan
        if os.path.isfile(f_out):
            try:
                om = read_outlier_map(f_out)
                pct_out = 100.0 * om.sum() / om.size
            except Exception:
                pass

        rows.append({
            "participant_id": sub,
            "diagnosis": r[dx_col],
            "pair_id": r[pair_col] if pair_col else np.nan,
            "n_volumes": len(rel_rms),
            "mean_abs_rms": float(abs_rms.mean()),
            "mean_rel_rms": float(rel_rms.mean()),
            "max_rel_rms": float(rel_rms.max()),
            "pct_outlier_slices": float(pct_out),
        })

    if missing:
        print(f"\n!! {len(missing)} subjects missing eddy output:")
        for s, p in missing[:10]:
            print(f"   {s}: {p}")
        print("   Fix these before continuing - do not silently proceed "
              "with a smaller cohort.")

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("FATAL: no motion data read. Check PREPROC and EDDY_STEM "
                 "paths at the top of this script.")

    # ---------------- exclusions ----------------
    df["fail_mean"] = df.mean_rel_rms > THRESH_MEAN_REL_RMS
    df["fail_max"] = df.max_rel_rms > THRESH_MAX_REL_RMS
    df["fail_outlier"] = df.pct_outlier_slices > THRESH_PCT_OUTLIER
    df["exclude"] = df.fail_mean | df.fail_max | df.fail_outlier

    if DROP_PARTNER and df.pair_id.notna().any():
        broken = set(df.loc[df.exclude, "pair_id"].dropna())
        df["exclude_final"] = df.exclude | df.pair_id.isin(broken)
    else:
        df["exclude_final"] = df.exclude

    # ---------------- report ----------------
    print("\n" + "=" * 68)
    print("MOTION SUMMARY BY GROUP")
    print("=" * 68)
    print(df.groupby("diagnosis")[
        ["mean_rel_rms", "max_rel_rms", "pct_outlier_slices"]
    ].agg(["mean", "std", "median", "max"]).round(3).to_string())

    groups = sorted(df.diagnosis.unique())
    if len(groups) == 2:
        a = df.loc[df.diagnosis == groups[0], "mean_rel_rms"].dropna()
        b = df.loc[df.diagnosis == groups[1], "mean_rel_rms"].dropna()
        t, pt = stats.ttest_ind(a, b, equal_var=False)
        u, pu = stats.mannwhitneyu(a, b)
        pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        d = (a.mean() - b.mean()) / pooled if pooled > 0 else np.nan

        print("\n" + "=" * 68)
        print(f"DOES MOTION DIFFER BY DIAGNOSIS?  "
              f"({groups[0]} n={len(a)} vs {groups[1]} n={len(b)})")
        print("=" * 68)
        print(f"  Welch t      t={t:+.3f}   p={pt:.4f}")
        print(f"  Mann-Whitney U={u:.0f}    p={pu:.4f}")
        print(f"  Cohen's d    d={d:+.3f}")
        if pt < 0.05:
            print("\n  >> Motion DIFFERS by group. It is a confound. It MUST")
            print("     go in the GLM as a covariate, and you must say so")
            print("     explicitly in the write-up.")
        else:
            print("\n  >> No significant group difference. Still include it")
            print("     as a covariate - absence of evidence at n=78 is not")
            print("     evidence of absence, and reviewers will ask.")

    print("\n" + "=" * 68)
    print("EXCLUSIONS")
    print("=" * 68)
    print(f"  thresholds: mean_rel_rms>{THRESH_MEAN_REL_RMS}  "
          f"max_rel_rms>{THRESH_MAX_REL_RMS}  "
          f"pct_outlier>{THRESH_PCT_OUTLIER}")
    bad = df[df.exclude]
    if bad.empty:
        print("  No subject exceeds any threshold. Cohort stays at "
              f"{len(df)}.")
    else:
        print(f"  {len(bad)} subject(s) fail directly:")
        for _, r in bad.iterrows():
            why = [n for n, f in (("mean", r.fail_mean),
                                  ("max", r.fail_max),
                                  ("outlier", r.fail_outlier)) if f]
            print(f"    {r.participant_id} ({r.diagnosis}) "
                  f"mean={r.mean_rel_rms:.2f} max={r.max_rel_rms:.2f} "
                  f"out={r.pct_outlier_slices:.1f}%  -> {','.join(why)}")
        n_extra = int(df.exclude_final.sum() - df.exclude.sum())
        if n_extra:
            print(f"  + {n_extra} partner(s) dropped to keep pairs intact.")
        print(f"  Cohort after exclusions: {int((~df.exclude_final).sum())}")

    df.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.4f")
    print(f"\nWrote {OUT_TSV}")
    print("Log the excluded subjects and the thresholds in "
          "cohort_decisions.md now, before running dtifit.")


if __name__ == "__main__":
    main()
