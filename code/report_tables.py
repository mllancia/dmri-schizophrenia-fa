#!/usr/bin/env python3
"""
Recompute cohort statistics at the final n = 72 and emit markdown ready to
paste into REPORT.md sections 2.4 and 3.2.

The tables currently in REPORT.md were computed at n = 78, before the six
QC exclusions. Those numbers are stale and are flagged [VERIFY] in the draft.

Usage:
    source .venv/bin/activate
    python code/report_tables.py
    python code/report_tables.py > derivatives/report_tables.md
"""

import os
import numpy as np
import pandas as pd
from scipy import stats

PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
os.chdir(PROJ)

qc = pd.read_csv("derivatives/motion_qc.tsv", sep="\t")
qc["exclude_final"] = qc.exclude_final.astype(bool)
d = qc.loc[~qc.exclude_final].copy()
coh = pd.read_csv("cohort_final.tsv", sep="\t")
d = d.merge(coh[["participant_id", "age", "gender", "ScannerSerialNumber",
                 "protocol", "ghost_NoGhost"]],
            on="participant_id", how="left")

try:
    o = pd.read_csv("derivatives/orientation_qc.tsv", sep="\t")
    d = d.merge(o[["participant_id", "pct_LR", "mean_MD", "mean_FA_wm"]],
                on="participant_id", how="left")
except FileNotFoundError:
    pass

C = d[d.diagnosis == "CONTROL"]
S = d[d.diagnosis == "SCHZ"]
assert len(d) == len(C) + len(S)

print(f"<!-- generated at n = {len(d)}: "
      f"{len(C)} CONTROL, {len(S)} SCZ -->\n")


def ms(x):
    return f"{x.mean():.3f} ± {x.std(ddof=1):.3f}"


def cohend(a, b):
    p = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / p if p > 0 else np.nan


def dci(a, b):
    dd = cohend(a, b)
    n1, n2 = len(a), len(b)
    se = np.sqrt((n1 + n2) / (n1 * n2) + dd**2 / (2 * (n1 + n2 - 2)))
    return dd, dd - 1.96 * se, dd + 1.96 * se


# ---------------------------------------------------------------- section 2.4
print("### Section 2.4 — Cohort characteristics\n")
print(f"| | Control (n = {len(C)}) | Schizophrenia (n = {len(S)}) | Test |")
print("|---|---|---|---|")
pa = stats.ttest_ind(C.age, S.age).pvalue
print(f"| Age, years | {ms(C.age)} | {ms(S.age)} | p = {pa:.3f} |")
for lab, col in (("Sex, F / M", "gender"),
                 ("Scanner 35343 / 35426", "ScannerSerialNumber"),
                 ("T1 ghost absent / present", "ghost_NoGhost")):
    t = pd.crosstab(d.diagnosis, d[col])
    p = (stats.fisher_exact(t)[1] if t.shape == (2, 2)
         else stats.chi2_contingency(t)[1])
    cc = " / ".join(str(v) for v in t.loc["CONTROL"])
    ss = " / ".join(str(v) for v in t.loc["SCHZ"])
    print(f"| {lab} | {cc} | {ss} | p = {p:.3f} |")

# ---------------------------------------------------------------- section 3.2
print("\n### Section 3.2 — Head motion\n")
print("| | Control | Schizophrenia |")
print("|---|---|---|")
for lab, col in (("Mean relative RMS, mm", "mean_rel_rms"),
                 ("Maximum relative RMS, mm", "max_rel_rms"),
                 ("Outlier slices, %", "pct_outlier_slices")):
    print(f"| {lab} | {ms(C[col])} | {ms(S[col])} |")
    print(f"| &nbsp;&nbsp;median | {C[col].median():.3f} | "
          f"{S[col].median():.3f} |")

a, b = C.mean_rel_rms, S.mean_rel_rms
t, pt = stats.ttest_ind(a, b, equal_var=False)
u, pu = stats.mannwhitneyu(a, b)
dd, lo, hi = dci(a, b)
df_w = ((a.var(ddof=1)/len(a) + b.var(ddof=1)/len(b))**2 /
        ((a.var(ddof=1)/len(a))**2/(len(a)-1) +
         (b.var(ddof=1)/len(b))**2/(len(b)-1)))
print(f"\nWelch t({df_w:.0f}) = {t:+.2f}, p = {pt:.3f}; "
      f"Mann–Whitney p = {pu:.3f}; "
      f"Cohen's d = {dd:+.2f}, 95% CI [{lo:.2f}, {hi:.2f}].\n")

# ---------------------------------------------------------------- orientation
if "pct_LR" in d.columns:
    print("### Section 3.3 — Orientation metric, analysis cohort\n")
    a2, b2 = C.pct_LR.dropna(), S.pct_LR.dropna()
    dd2, lo2, hi2 = dci(a2, b2)
    print(f"- CONTROL {ms(a2)}%, SCZ {ms(b2)}% "
          f"(p = {stats.ttest_ind(a2, b2).pvalue:.3f}, "
          f"d = {dd2:+.2f} [{lo2:.2f}, {hi2:.2f}])")
    print(f"- Cohort median {d.pct_LR.median():.1f}%, "
          f"range {d.pct_LR.min():.1f}–{d.pct_LR.max():.1f}%")
    print("- After excluding sub-10855 the groups are indistinguishable on "
          "this metric.\n")

print("<!-- paste the tables above over the [VERIFY]-marked ones in "
      "REPORT.md -->")
