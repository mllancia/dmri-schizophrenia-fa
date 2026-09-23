#!/usr/bin/env python3
"""
Phase 5 - Colour-FA orientation check.

THE POINT
---------
This is the ONLY reliable defence against the silent rotated-bvec bug. If
dtifit was given the wrong gradient directions, FA still looks perfectly
normal - bright corpus callosum, dark ventricles, plausible values. Only the
ORIENTATIONS are wrong, and the FA differences that follow look like real
findings.

WHAT YOU ARE LOOKING FOR
------------------------
Colour encodes the direction of V1, the principal eigenvector:
    RED   = left-right          -> corpus callosum (crosses the midline)
    BLUE  = superior-inferior   -> corticospinal tract (runs up the brainstem)
    GREEN = anterior-posterior  -> superior longitudinal fasciculus

On a mid-axial slice you should see a clear RED band front and back on the
midline (genu and splenium of the corpus callosum), GREEN running front-to-back
either side of it, and BLUE where descending fibres pass through.

If the colours are swapped, scrambled, or vary randomly between subjects,
STOP. Do not proceed to ROI analysis.

Usage:
    source ~/VeronaSemester/DmriMiniCourse/.venv/bin/activate
    python code/phase5_colorfa_check.py            # montage of all subjects
    python code/phase5_colorfa_check.py sub-10159  # detail for one subject
"""

import os
import sys
import glob
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = os.path.expanduser("~/VeronaSemester/DmriMiniCourse")
DTIDIR = os.path.join(PROJ, "derivatives", "dti")
QCDIR = os.path.join(PROJ, "derivatives", "qc")
os.makedirs(QCDIR, exist_ok=True)


def colour_fa(sub_dir):
    """Return an (x, y, z, 3) RGB array = |V1| * FA, clipped to [0, 1]."""
    fa = nib.load(os.path.join(sub_dir, "dti_FA.nii.gz")).get_fdata()
    v1 = nib.load(os.path.join(sub_dir, "dti_V1.nii.gz")).get_fdata()
    rgb = np.abs(v1) * fa[..., None]
    return np.clip(np.nan_to_num(rgb), 0, 1)


def axial(rgb, frac):
    """One axial slice, oriented for display (radiological-ish)."""
    z = int(rgb.shape[2] * frac)
    return np.transpose(rgb[:, :, z, :], (1, 0, 2))[::-1]


def montage():
    subs = sorted(d for d in glob.glob(os.path.join(DTIDIR, "sub-*"))
                  if os.path.isfile(os.path.join(d, "dti_V1.nii.gz")))
    if not subs:
        sys.exit(f"No dti_V1 maps found in {DTIDIR}. Run phase5_dtifit.sh.")

    n = len(subs)
    ncol = 9
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 1.7, nrow * 1.9),
                             facecolor="black")
    axes = np.atleast_2d(axes).ravel()

    for ax in axes:
        ax.axis("off")

    for ax, d in zip(axes, subs):
        try:
            ax.imshow(axial(colour_fa(d), 0.50), interpolation="nearest")
            ax.set_title(os.path.basename(d).replace("sub-", ""),
                         color="white", fontsize=6, pad=2)
        except Exception as e:
            ax.set_title(f"ERR {os.path.basename(d)}", color="red", fontsize=5)
            print(f"  failed {d}: {e}")

    out = os.path.join(QCDIR, "colorfa_montage.png")
    fig.suptitle(f"Colour-FA orientation check  |  {n} subjects  |  "
                 "CC should be RED, CST BLUE, SLF GREEN",
                 color="white", fontsize=11, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    fig.savefig(out, dpi=140, facecolor="black")
    print(f"Wrote {out}")
    print(f"\nOpen it and scan every panel. You are checking two things:")
    print("  1. Each brain individually shows red CC / green SLF / blue CST.")
    print("  2. All 78 look the SAME. One odd subject means one broken fit,")
    print("     which is arguably worse than all of them being broken,")
    print("     because it will not be obvious downstream.")


def detail(sub):
    if not sub.startswith("sub-"):
        sub = "sub-" + sub
    d = os.path.join(DTIDIR, sub)
    if not os.path.isdir(d):
        sys.exit(f"No such subject directory: {d}")

    rgb = colour_fa(d)
    fig, axes = plt.subplots(1, 3, figsize=(13, 5), facecolor="black")
    for ax, frac, lab in zip(axes, (0.40, 0.50, 0.60),
                             ("lower", "mid", "upper")):
        ax.imshow(axial(rgb, frac), interpolation="nearest")
        ax.set_title(f"axial {lab} ({frac:.0%} of z)", color="white",
                     fontsize=10)
        ax.axis("off")

    fig.suptitle(f"{sub}  -  RED=L-R (corpus callosum)   "
                 "GREEN=A-P (SLF)   BLUE=S-I (corticospinal)",
                 color="white", fontsize=11)
    out = os.path.join(QCDIR, f"colorfa_{sub}.png")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out, dpi=150, facecolor="black")
    print(f"Wrote {out}")

    fa = nib.load(os.path.join(d, "dti_FA.nii.gz")).get_fdata()
    inb = fa[fa > 0.05]
    print(f"\n{sub} FA: max={fa.max():.3f}  "
          f"mean(in-brain)={inb.mean():.3f}  "
          f"frac>0.2={(inb > 0.2).mean():.2f}")
    print("Expect max just under 1.0, in-brain mean roughly 0.25-0.40.")
    if fa.max() > 1.0:
        print("!! FA exceeds 1.0 - this is impossible. Investigate before "
              "going further.")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        detail(sys.argv[1])
    else:
        montage()
