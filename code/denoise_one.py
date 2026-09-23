import sys, time, os
import numpy as np, nibabel as nib
from dipy.denoise.localpca import mppca
from dipy.denoise.gibbs import gibbs_removal


def main():
    sub = sys.argv[1]
    src = f'raw/{sub}/dwi/{sub}_dwi.nii.gz'
    out = f'derivatives/preproc/{sub}'
    os.makedirs(out, exist_ok=True)
    os.makedirs('qc', exist_ok=True)

    if os.path.exists(f'{out}/{sub}_degibbs.nii.gz'):
        print(f'{sub}  already done, skipping', flush=True)
        return

    img = nib.load(src)
    data = img.get_fdata(dtype=np.float32)
    print(f'{sub}  shape {data.shape}', flush=True)
    assert data.ndim == 4 and data.shape[3] == 65, 'unexpected shape'

    t = time.time()
    den, sigma = mppca(data, patch_radius=2, return_sigma=True)
    print(f'  mppca         {time.time()-t:7.1f} s', flush=True)
    nib.save(nib.Nifti1Image(sigma.astype(np.float32), img.affine),
             f'qc/{sub}_noise.nii.gz')

    t = time.time()
    deg = gibbs_removal(den, slice_axis=2, num_processes=4)
    print(f'  gibbs_removal {time.time()-t:7.1f} s', flush=True)
    nib.save(nib.Nifti1Image(deg.astype(np.float32), img.affine),
             f'{out}/{sub}_degibbs.nii.gz')
    print('  done', flush=True)


if __name__ == '__main__':
    main()
