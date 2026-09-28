# Diffusion MRI: White Matter Microstructure in Schizophrenia

Course project (Diffusion MRI mini-course, University of Verona) testing white
matter fractional anisotropy differences between schizophrenia and control
participants in OpenNeuro dataset ds000030 (UCLA CNP LA5c).

Analysis cohort: 72 participants (36 schizophrenia, 36 control), matched
exactly on acquisition protocol and sex and approximately on age.

## Contents

| Path | Contents |
|---|---|
| `REPORT.md` | Full write-up: methods, QC, results, discussion |
| `cohort_decisions.md` | Every methodological decision with its evidence |
| `code/` | Preprocessing, QC, tensor fitting, registration, analysis |
| `main.ipynb` | Cohort construction and protocol-confound analysis |
| `validate_fit.py` | From-scratch tensor fit, used as a positive control |
| `cohort_final.tsv` | Cohort definition with acquisition metadata |
| `matched_pairs_final.tsv` | Matched pairings (analysis is unpaired by design) |
| `derivatives/*.tsv` | QC metrics and regional statistics |
| `docs/ds000030_dataset_README.md` | Original dataset README, kept as source evidence |

Imaging data are not included. They are public and can be retrieved from
OpenNeuro with `fetch_dwi.sh`.

## Pipeline

MP-PCA denoising → Gibbs removal → BET → eddy (`--repol`) → dtifit →
nonlinear registration to FMRIB58_FA → JHU atlas ROI analysis → TBSS with
permutation testing

Run order: `run_batch.sh` → `phase4f_motion_qc.py` → `phase5_dtifit.sh` →
`phase5_colorfa_check.py` → `phase6a_register.sh` → `phase6b_roi_extract.py` →
`phase6c_interrogate.py` → `phase8_tbss.sh` → `phase8b_tbss_report.py`

All analysis subject lists derive from `derivatives/motion_qc.tsv` filtered on
`exclude_final`, never from directory listings: `raw/` holds ~110 participants
and `derivatives/dti/` holds 74, while the analysis cohort is 72.

## Requirements

FSL, plus the Python packages in `requirements.txt`.

## Data citation

Poldrack RA, Congdon E, Triplett W, et al. A phenome-wide examination of neural
and cognitive function. Scientific Data 2016;3:160110.
doi:10.1038/sdata.2016.110
