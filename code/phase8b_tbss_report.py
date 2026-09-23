#!/usr/bin/env python3
"""
Phase 8b - Summarise TBSS results for the report.

Produces:
  - extent of significant skeleton for each contrast (1-p > 0.95)
  - which JHU tracts the significant voxels fall in, as % of each tract's
    skeleton - this is what tells you whether the result is DIFFUSE (many
    tracts partially significant, consistent with the global shift in §4.2)
    or FOCAL (a few tracts heavily significant)
  - a figure: significant skeleton over the mean FA, several axial slices

Usage:
    source .venv/bin/activate
    python code/phase8b_tbss_report.py
"""

import os
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
os.chdir(PROJ)
FSLDIR = os.environ.get("FSLDIR", os.path.expanduser("~/fsl"))
S = "derivatives/tbss/stats"


def load(p):
    return nib.load(p).get_fdata()


mean = load(f"{S}/mean_FA.nii.gz")
skel = load(f"{S}/mean_FA_skeleton_mask.nii.gz") > 0
atlas = load(f"{FSLDIR}/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz").astype(int)

names = {}
try:
    import xml.etree.ElementTree as ET
    for lab in ET.parse(f"{FSLDIR}/data/atlases/JHU-labels.xml").getroot().iter("label"):
        names[int(lab.get("index")) + 1] = lab.text.strip()
except Exception:
    pass

print(f"Skeleton voxels: {skel.sum()}\n")

for t, lab in ((1, "CONTROL > SCZ"), (2, "SCZ > CONTROL")):
    f = f"{S}/tbss_FA_tfce_corrp_tstat{t}.nii.gz"
    if not os.path.isfile(f):
        print(f"missing {f}")
        continue
    p = load(f)
    sig = (p > 0.95) & skel
    print("=" * 70)
    print(f"Contrast {t}: {lab}")
    print("=" * 70)
    print(f"  max 1-p = {p[skel].max():.4f}")
    print(f"  significant voxels: {sig.sum()} ({100*sig.sum()/skel.sum():.1f}% "
          f"of skeleton)")
    if sig.sum() == 0:
        print("  No significant voxels.\n")
        continue

    rows = []
    for l in sorted(set(np.unique(atlas[skel])) - {0}):
        tract_sk = skel & (atlas == l)
        n = int(tract_sk.sum())
        if n < 20:
            continue
        k = int((sig & (atlas == l)).sum())
        rows.append((names.get(l, f"label_{l}"), k, n, 100 * k / n))
    rows.sort(key=lambda r: -r[3])
    unl = int((sig & (atlas == 0)).sum())
    print(f"\n  {'JHU tract':<58}{'sig':>6}{'skel':>6}{'%':>7}")
    for nm, k, n, pc in rows[:20]:
        if k:
            print(f"  {nm[:57]:<58}{k:>6}{n:>6}{pc:>6.1f}%")
    print(f"  (outside labelled JHU regions: {unl} voxels)")
    ntr = sum(1 for r in rows if r[1] > 0)
    print(f"\n  Tracts with any significant voxels: {ntr} of {len(rows)}")
    print("  Many tracts partially significant => diffuse, consistent with the")
    print("  global shift. Few tracts heavily significant => focal.\n")

    fig, ax = plt.subplots(2, 4, figsize=(15, 8), facecolor="black")
    zs = np.linspace(0.30, 0.72, 8) * mean.shape[2]
    for a, z in zip(ax.ravel(), zs.astype(int)):
        a.imshow(np.rot90(mean[:, :, z]), cmap="gray", vmin=0, vmax=0.8)
        sk = np.rot90(skel[:, :, z]).astype(float)
        sg = np.rot90(sig[:, :, z]).astype(float)
        a.imshow(np.ma.masked_where(sk == 0, sk), cmap="summer", alpha=0.6)
        a.imshow(np.ma.masked_where(sg == 0, sg), cmap="autumn", alpha=0.95)
        a.set_title(f"z = {z}", color="w", fontsize=9)
        a.axis("off")
    fig.suptitle(f"TBSS {lab}: green = skeleton, red/yellow = TFCE-corrected "
                 f"1-p > 0.95  ({sig.sum()} voxels)", color="w", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = f"derivatives/figures_tbss_tstat{t}.png"
    plt.savefig(out, dpi=140, facecolor="black")
    print(f"  wrote {out}\n")
