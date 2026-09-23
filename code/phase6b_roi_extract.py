#!/usr/bin/env python3
"""
Phase 6b - JHU atlas ROI extraction and the primary hypothesis test.

WHAT THIS DOES
--------------
1. Applies the JHU ICBM-DTI-81 atlas to every subject's standard-space FA, MD,
   AD and RD maps, excluding voxels with FA > 1 (degenerate OLS fits).
2. Tests the PRE-SPECIFIED primary hypothesis: reduced corpus callosum FA in
   schizophrenia, adjusted for age, sex and motion.
3. Runs the exploratory regional analysis across all remaining JHU regions,
   with FDR correction, clearly labelled exploratory.

Effect sizes and confidence intervals are reported alongside p-values
throughout. At n = 72 the p-value alone is close to uninformative.

Usage:
    source .venv/bin/activate
    pip install statsmodels        # if not already present
    python code/phase6b_roi_extract.py
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import nibabel as nib

PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
os.chdir(PROJ)
FSLDIR = os.environ.get("FSLDIR", os.path.expanduser("~/fsl"))

ATLAS = f"{FSLDIR}/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz"
LABELS_XML = f"{FSLDIR}/data/atlases/JHU-labels.xml"
STDDIR = "derivatives/std"
OUT_LONG = "derivatives/roi_values.tsv"
OUT_STATS = "derivatives/roi_stats.tsv"

# Pre-specified primary ROI: corpus callosum (genu + body + splenium).
CC_LABELS = {3: "Genu of corpus callosum",
             4: "Body of corpus callosum",
             5: "Splenium of corpus callosum"}

if not os.path.isfile(ATLAS):
    sys.exit(f"FATAL: JHU atlas not found at {ATLAS}\n"
             "Check FSLDIR, or locate it with:\n"
             "  find $FSLDIR/data/atlases -name 'JHU-ICBM-labels*'")


def load_labels():
    """Parse region names from FSL's atlas XML; fall back to numeric names."""
    try:
        import xml.etree.ElementTree as ET
        root = ET.parse(LABELS_XML).getroot()
        out = {}
        for lab in root.iter("label"):
            out[int(lab.get("index")) + 1] = lab.text.strip()
        return out
    except Exception as e:
        print(f"  (could not parse {LABELS_XML}: {e}; using numeric names)")
        return {}


def main():
    qc = pd.read_csv("derivatives/motion_qc.tsv", sep="\t")
    qc["exclude_final"] = qc.exclude_final.astype(bool)
    keep = qc.loc[~qc.exclude_final].copy()
    coh = pd.read_csv("cohort_final.tsv", sep="\t")
    keep = keep.merge(coh[["participant_id", "age", "gender"]],
                      on="participant_id", how="left")
    print(f"Analysis cohort: {len(keep)} subjects "
          f"({(keep.diagnosis=='CONTROL').sum()} CONTROL, "
          f"{(keep.diagnosis=='SCHZ').sum()} SCZ)")

    atlas = nib.load(ATLAS).get_fdata().astype(int)
    names = load_labels()
    present = sorted(set(np.unique(atlas)) - {0})
    print(f"Atlas regions present: {len(present)}")

    rows, missing = [], []
    for _, r in keep.iterrows():
        sub = r.participant_id
        d = f"{STDDIR}/{sub}"
        if not os.path.isfile(f"{d}/FA_std.nii.gz"):
            missing.append(sub)
            continue
        maps = {}
        for m in ("FA", "MD", "AD", "RD"):
            f = f"{d}/{m}_std.nii.gz"
            maps[m] = nib.load(f).get_fdata() if os.path.isfile(f) else None
        if maps["FA"] is None:
            missing.append(sub)
            continue
        if maps["FA"].shape != atlas.shape:
            sys.exit(f"FATAL: {sub} FA_std shape {maps['FA'].shape} != atlas "
                     f"{atlas.shape}. Registration target mismatch.")

        # THE RULE: degenerate fits are excluded from every average.
        valid = (maps["FA"] > 0) & (maps["FA"] <= 1)

        for lab in present:
            m = (atlas == lab) & valid
            n = int(m.sum())
            if n < 10:
                continue
            rec = {"participant_id": sub, "diagnosis": r.diagnosis,
                   "age": r.age, "gender": r.gender,
                   "motion": r.mean_rel_rms,
                   "label": lab, "region": names.get(lab, f"label_{lab}"),
                   "n_voxels": n}
            for k, v in maps.items():
                rec[k] = float(v[m].mean()) if v is not None else np.nan
            rows.append(rec)

    if missing:
        print(f"\n!! {len(missing)} subjects missing standard-space maps: "
              f"{missing[:6]}{'...' if len(missing) > 6 else ''}")
        print("   Rerun code/phase6a_register.sh before trusting these results.")

    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("FATAL: no ROI values extracted.")
    df.to_csv(OUT_LONG, sep="\t", index=False, float_format="%.6f")
    print(f"Wrote {OUT_LONG}  ({df.participant_id.nunique()} subjects x "
          f"{df.label.nunique()} regions)")

    try:
        import statsmodels.formula.api as smf
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        sys.exit("\nInstall statsmodels to run the models:\n"
                 "  pip install statsmodels")

    def model(sub_df, measure):
        """Adjusted group difference: measure ~ diagnosis + age + sex + motion."""
        t = sub_df.dropna(subset=[measure, "age", "gender", "motion"]).copy()
        t["dx"] = (t.diagnosis == "SCHZ").astype(int)
        t["sex"] = (t.gender.astype(str).str.upper().str[0] == "M").astype(int)
        fit = smf.ols(f"{measure} ~ dx + age + sex + motion", data=t).fit()
        b = fit.params["dx"]
        ci = fit.conf_int().loc["dx"]
        a = t.loc[t.dx == 0, measure]
        c = t.loc[t.dx == 1, measure]
        pooled = np.sqrt((a.var(ddof=1) + c.var(ddof=1)) / 2)
        return {"n": len(t), "control_mean": a.mean(), "scz_mean": c.mean(),
                "beta_dx": b, "ci_lo": ci[0], "ci_hi": ci[1],
                "p": fit.pvalues["dx"],
                "cohen_d": (c.mean() - a.mean()) / pooled if pooled else np.nan}

    # ------------------------------------------------- PRIMARY
    print("\n" + "=" * 70)
    print("PRIMARY HYPOTHESIS (pre-specified): corpus callosum FA")
    print("=" * 70)
    cc = df[df.label.isin(CC_LABELS)]
    if cc.empty:
        print("  !! No corpus callosum labels found. Check atlas indices "
              "against your FSL version:")
        print(df[["label", "region"]].drop_duplicates().head(12).to_string(index=False))
    else:
        agg = (cc.groupby(["participant_id", "diagnosis", "age", "gender",
                           "motion"], as_index=False)
                 .apply(lambda g: pd.Series({
                     "FA": np.average(g.FA, weights=g.n_voxels)}))
                 .reset_index(drop=True))
        res = model(agg, "FA")
        print(f"  n = {res['n']}")
        print(f"  CONTROL mean FA  {res['control_mean']:.4f}")
        print(f"  SCZ     mean FA  {res['scz_mean']:.4f}")
        print(f"  adjusted difference  {res['beta_dx']:+.5f}  "
              f"95% CI [{res['ci_lo']:+.5f}, {res['ci_hi']:+.5f}]")
        print(f"  Cohen's d  {res['cohen_d']:+.3f}")
        print(f"  p = {res['p']:.4f}   {'<-- significant' if res['p']<0.05 else ''}")
        print("\n  Report the effect size and CI first. At ~40% power a "
              "non-significant\n  result constrains the effect size; it does "
              "not demonstrate absence.")
        print("\n  Subregions:")
        for lab, nm in CC_LABELS.items():
            s = df[df.label == lab]
            if not s.empty:
                r = model(s, "FA")
                print(f"    {nm:<32} d={r['cohen_d']:+.3f}  p={r['p']:.4f}")

    # ------------------------------------------------- EXPLORATORY
    print("\n" + "=" * 70)
    print("EXPLORATORY: all JHU regions, FA")
    print("=" * 70)
    out = []
    for lab, g in df.groupby("label"):
        if g.participant_id.nunique() < 20:
            continue
        r = model(g, "FA")
        r.update({"label": lab, "region": g.region.iloc[0]})
        out.append(r)
    st = pd.DataFrame(out)
    st["p_fdr"] = multipletests(st.p, method="fdr_bh")[1]
    st = st.sort_values("p")
    st.to_csv(OUT_STATS, sep="\t", index=False, float_format="%.6f")

    print(st[["region", "control_mean", "scz_mean", "cohen_d", "p", "p_fdr"]]
          .head(15).to_string(index=False,
                              float_format=lambda x: f"{x:.4f}"))
    n_sig = int((st.p_fdr < 0.05).sum())
    print(f"\nRegions surviving FDR q < 0.05: {n_sig} of {len(st)}")
    print(f"Wrote {OUT_STATS}")
    print("\nEVERYTHING IN THIS SECTION IS EXPLORATORY. Label it as such in "
          "the report,\nand do not present any region here as a confirmed "
          "finding.")


if __name__ == "__main__":
    main()
