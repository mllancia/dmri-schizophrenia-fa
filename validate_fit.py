"""
Single-voxel tensor fit, used to cross-check FSL's dtifit on this study's data.

Lifted unchanged from the Brusini course project (fit_dti.py): design_matrix,
fit_dti, and tensor_metrics are verbatim. Everything ROI-specific has been
dropped. This is a verification tool, not a pipeline -- dtifit produces the
whole-brain maps.

Usage:
    python validate_fit.py <dwi> <bval> <rotated_bvec> <i> <j> <k>

Example (comparing against dtifit output at the same voxel):
    python validate_fit.py eddy_out.nii.gz sub-10159_dwi.bval \
        eddy_out.eddy_rotated_bvecs 48 55 40
"""

import sys

import numpy as np
from dipy.io import read_bvals_bvecs
from dipy.io.image import load_nifti


# --------------------------------------------------------------------------
# Verbatim from fit_dti.py
# --------------------------------------------------------------------------
def design_matrix(bvals, bvecs):
    """Build the (N, 6) matrix A from b-values and gradient directions.

    Each row is  A_i = -b_i * [ux^2, uy^2, uz^2, 2*ux*uy, 2*ux*uz, 2*uy*uz],
    matching the coefficient order d = [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz].
    b = 0 rows become all-zeros and contribute nothing to the fit.
    """
    ux, uy, uz = bvecs[:, 0], bvecs[:, 1], bvecs[:, 2]
    A = np.column_stack(
        [ux**2, uy**2, uz**2, 2 * ux * uy, 2 * ux * uz, 2 * uy * uz]
    )
    return -bvals[:, None] * A


def fit_dti(signal, bvals, bvecs):
    """Estimate the diffusion tensor for a single voxel. Returns (3, 3) D."""
    signal = np.asarray(signal, dtype=float)

    b0_mask = bvals == 0
    eps = 1e-8
    s0 = max(signal[b0_mask].mean(), eps)

    y = np.log(np.clip(signal, eps, None) / s0)

    A = design_matrix(bvals, bvecs)
    d = np.linalg.pinv(A) @ y

    Dxx, Dyy, Dzz, Dxy, Dxz, Dyz = d
    return np.array(
        [
            [Dxx, Dxy, Dxz],
            [Dxy, Dyy, Dyz],
            [Dxz, Dyz, Dzz],
        ]
    )


def tensor_metrics(D):
    """Return (FA, MD, AD, RD) from a 3x3 diffusion tensor."""
    evals = np.linalg.eigh(D)[0][::-1]
    l1, l2, l3 = evals

    ad = l1
    rd = (l2 + l3) / 2
    md = (l1 + l2 + l3) / 3

    num = (l1 - md) ** 2 + (l2 - md) ** 2 + (l3 - md) ** 2
    den = l1**2 + l2**2 + l3**2
    fa = np.sqrt(1.5) * np.sqrt(num / den) if den > 0 else 0.0
    fa = float(np.clip(fa, 0.0, 1.0))

    return fa, md, ad, rd


# --------------------------------------------------------------------------
# New: single-voxel entry point
# --------------------------------------------------------------------------
def fit_voxel(dwi_file, bval_file, bvec_file, i, j, k):
    data, _ = load_nifti(dwi_file)
    bvals, bvecs = read_bvals_bvecs(bval_file, bvec_file)

    if data.shape[3] != bvals.size:
        raise SystemExit(
            f"Volume/gradient mismatch: image has {data.shape[3]} volumes, "
            f"bval file has {bvals.size} entries. Stop and investigate."
        )

    signal = data[i, j, k, :]
    D = fit_dti(signal, bvals, bvecs)
    fa, md, ad, rd = tensor_metrics(D)

    print(f"voxel ({i}, {j}, {k})  of image {data.shape}")
    print(f"  S0 (mean b0) = {signal[bvals == 0].mean():.1f}")
    print(f"  FA = {fa:.4f}")
    print(f"  MD = {md:.6e}")
    print(f"  AD = {ad:.6e}")
    print(f"  RD = {rd:.6e}")
    print()
    print("Compare against dtifit with:")
    print(f"  fslmeants -i <prefix>_FA -c {i} {j} {k}")
    print(f"  fslmeants -i <prefix>_MD -c {i} {j} {k}")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    dwi, bval, bvec = sys.argv[1:4]
    i, j, k = (int(v) for v in sys.argv[4:7])
    fit_voxel(dwi, bval, bvec, i, j, k)
