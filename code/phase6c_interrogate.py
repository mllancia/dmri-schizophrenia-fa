#!/usr/bin/env python3
"""
Phase 6c - Interrogate the exploratory finding before writing it up.

Nine regions survive FDR and every effect points the same way. That pattern
has two very different explanations:

  (A) a genuinely widespread white matter difference, or
  (B) a systematic quality/SNR difference between groups that depresses FA
      everywhere and has nothing to do with biology.

These are indistinguishable from the regional table alone. This script runs
the checks that separate them, plus the AD/RD characterisation and a forest
plot for the report.

Usage:
    source .venv/bin/activate
    python code/phase6c_interrogate.py
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
os.chdir(PROJ)

df = pd.read_csv("derivatives/roi_values.tsv", sep="\t")
st = pd.read_csv("derivatives/roi_stats.tsv", sep="\t")
df["dx"] = (df.diagnosis == "SCHZ").astype(int)
df["sex"] = (df.gender.astype(str).str.upper().str[0] == "M").astype(int)


def adj(t, measure):
    """Adjusted group effect with Cohen's d and its CI."""
    t = t.dropna(subset=[measure, "age", "sex", "motion"])
    f = smf.ols(f"{measure} ~ dx + age + sex + motion", data=t).fit()
    a = t.loc[t.dx == 0, measure]
    c = t.loc[t.dx == 1, measure]
    p = np.sqrt((a.var(ddof=1) + c.var(ddof=1)) / 2)
    d = (c.mean() - a.mean()) / p if p else np.nan
    n1, n2 = len(a), len(c)
    se = np.sqrt((n1+n2)/(n1*n2) + d**2/(2*(n1+n2-2)))
    return d, d-1.96*se, d+1.96*se, f.pvalues["dx"], f.params["dx"]


# ============================================================ 1. GLOBAL
print("=" * 72)
print("1. IS THIS A GLOBAL SHIFT OR REGIONAL?")
print("=" * 72)
g = (df.groupby(["participant_id","diagnosis","dx","age","sex","motion"],
                as_index=False)
       .apply(lambda x: pd.Series({
           "FA": np.average(x.FA, weights=x.n_voxels),
           "MD": np.average(x.MD, weights=x.n_voxels)}))
       .reset_index(drop=True))
d, lo, hi, p, b = adj(g, "FA")
print(f"  Whole-atlas mean FA: CONTROL {g.loc[g.dx==0,'FA'].mean():.4f}  "
      f"SCZ {g.loc[g.dx==1,'FA'].mean():.4f}")
print(f"  adjusted diff {b:+.5f}   d = {d:+.3f} [{lo:.2f}, {hi:.2f}]   "
      f"p = {p:.4f}")
print("\n  If this global effect is significant and comparable in size to the")
print("  strongest regional effects, the finding is a WHOLE-BRAIN shift, not")
print("  a set of specific tracts. Say so plainly rather than listing regions")
print("  as though they were selectively affected.")

print("\n  Direction consistency across all regions:")
neg = int((st.cohen_d < 0).sum())
sp = stats.binomtest(neg, len(st), 0.5).pvalue
print(f"    {neg} of {len(st)} regions show lower FA in SCZ "
      f"(sign test p = {sp:.2g})")
print("    CAVEAT: regions are spatially contiguous and highly correlated,")
print("    so this p-value is anticonservative. Report it descriptively.")

# ==================================================== 2. ARTIFACT OR BIOLOGY
print("\n" + "=" * 72)
print("2. COULD A QUALITY DIFFERENCE EXPLAIN IT?")
print("=" * 72)
sub = df.drop_duplicates("participant_id")[
    ["participant_id","dx","motion","age","sex"]].merge(
    g[["participant_id","FA"]], on="participant_id")
try:
    o = pd.read_csv("derivatives/orientation_qc.tsv", sep="\t")
    sub = sub.merge(o[["participant_id","mean_MD","pct_LR"]],
                    on="participant_id", how="left")
except FileNotFoundError:
    pass
q = pd.read_csv("derivatives/motion_qc.tsv", sep="\t")
sub = sub.merge(q[["participant_id","pct_outlier_slices","max_rel_rms"]],
                on="participant_id", how="left")

for col, label in (("motion","mean relative displacement"),
                   ("pct_outlier_slices","% outlier slices"),
                   ("max_rel_rms","max displacement")):
    if col not in sub: continue
    a, c = sub.loc[sub.dx==0,col], sub.loc[sub.dx==1,col]
    r = np.corrcoef(sub[col].fillna(sub[col].mean()), sub.FA)[0,1]
    print(f"  {label:<32} group p = {stats.ttest_ind(a,c).pvalue:.3f}   "
          f"corr with global FA = {r:+.3f}")

print("\n  A quality confound would show BOTH a group difference in the")
print("  quality metric AND a strong correlation with global FA. If neither")
print("  holds, the artifact explanation is weak and you can argue for")
print("  interpretation (A). If motion correlates strongly with FA, say so")
print("  even though it is already a covariate.")

# ==================================================== 3. OUTLIER-DRIVEN?
print("\n" + "=" * 72)
print("3. IS THE TOP RESULT DRIVEN BY A FEW SUBJECTS?")
print("=" * 72)
top = st.sort_values("p").iloc[0]
t = df[df.region == top.region].copy()
d0 = adj(t, "FA")[0]
loo = []
for s in t.participant_id:
    dd = adj(t[t.participant_id != s], "FA")[0]
    loo.append((s, dd))
loo.sort(key=lambda x: abs(x[1] - d0), reverse=True)
print(f"  {top.region[:60]}")
print(f"  full-sample d = {d0:+.3f}")
print("  most influential subjects (leave-one-out d):")
for s, dd in loo[:5]:
    print(f"    without {s}: d = {dd:+.3f}  (shift {dd-d0:+.3f})")
print("\n  A stable effect barely moves. If removing one subject changes d by")
print("  more than ~0.15, the finding rests on that subject and must be")
print("  reported with that caveat.")

# ==================================================== 4. AD / RD
print("\n" + "=" * 72)
print("4. AD AND RD IN THE FDR-SIGNIFICANT REGIONS")
print("=" * 72)
sig = st[st.p_fdr < 0.05].region.tolist()
print(f"{'region':<46}{'FA d':>7}{'AD d':>7}{'RD d':>7}  pattern")
print("-" * 82)
for rg in sig:
    t = df[df.region == rg]
    ds = {m: adj(t, m)[0] for m in ("FA","AD","RD")}
    if ds["RD"] > 0.2 and abs(ds["AD"]) < 0.2:
        pat = "RD up only"
    elif ds["AD"] < -0.2 and abs(ds["RD"]) < 0.2:
        pat = "AD down only"
    elif ds["RD"] > 0.2 and ds["AD"] < -0.2:
        pat = "both"
    else:
        pat = "unclear"
    print(f"{rg[:45]:<46}{ds['FA']:>7.2f}{ds['AD']:>7.2f}{ds['RD']:>7.2f}  {pat}")
print("\n  CAUTION: AD and RD are no more biologically specific than FA.")
print("  Describe the pattern; do not claim demyelination or axonal loss.")

# ==================================================== 5. FIGURE
print("\n" + "=" * 72)
print("5. FOREST PLOT")
print("=" * 72)
rows = []
for rg in st.sort_values("cohen_d").region.head(20):
    t = df[df.region == rg]
    d, lo, hi, p, _ = adj(t, "FA")
    rows.append((rg, d, lo, hi, st.loc[st.region==rg,"p_fdr"].iloc[0]))
rows.reverse()
fig, ax = plt.subplots(figsize=(10, 8))
for i, (rg, d, lo, hi, pf) in enumerate(rows):
    c = "#C2410C" if pf < 0.05 else "#94A3B8"
    ax.plot([lo, hi], [i, i], color=c, lw=2)
    ax.plot(d, i, "o", color=c, ms=6)
ax.axvline(0, color="black", lw=1, ls="--")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0][:48] for r in rows], fontsize=8)
ax.set_xlabel("Cohen's d (negative = lower FA in schizophrenia)")
ax.set_title("Regional FA differences, JHU atlas (exploratory)\n"
             "coloured = FDR q < 0.05", fontsize=11)
plt.tight_layout()
os.makedirs("derivatives/qc", exist_ok=True)
plt.savefig("derivatives/figures_forest_FA.png", dpi=150)
print("  wrote derivatives/figures_forest_FA.png")
