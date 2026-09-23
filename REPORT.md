# White Matter Microstructure in Schizophrenia: A Diffusion Tensor Imaging Study of the UCLA Consortium for Neuropsychiatric Phenomics Cohort

**Frederick Lancia**
Erasmus Mundus Maths DISC — University of Verona
Diffusion MRI Mini-Course Project

---

> **STATUS OF THIS DRAFT**
> Sections 1–3 and 6–7 are complete. Sections 4 (group results) and 5
> (discussion) contain `[TODO]` blocks that can only be filled once Phases 6–8
> have run. Everything marked `[VERIFY]` is a number or fact that should be
> confirmed against the pipeline before submission.
> Delete this box before submitting.

---

## Abstract

Reduced fractional anisotropy (FA) in white matter is among the most replicated
neuroimaging findings in schizophrenia, and is commonly interpreted as evidence
for the disconnection hypothesis. The disconnectino hypothesis states that psychiatric and neurological disorders stem from a breakdown in the functional integration and communication between distributed brain regions rather than damage to isolated centers. We tested for group differences in white
matter FA between individuals with schizophrenia and healthy controls using
diffusion-weighted imaging from the UCLA Consortium for Neuropsychiatric
Phenomics LA5c dataset (OpenNeuro ds000030).

From 272 available participants we constructed a cohort matched exactly on
acquisition protocol and sex and approximately on age. This matching was
necessary because acquisition protocol was strongly confounded with diagnosis
in the full sample (χ², p = 0.0001), an imbalance that would have rendered any
FA finding uninterpretable. After preprocessing and a two-stage quality control
procedure, 72 participants (36 schizophrenia, 36 control) entered analysis.

`[TODO — one sentence stating the primary corpus callosum result with effect
size and CI]`
`[TODO — one sentence on exploratory ROI and TBSS findings]`

We additionally report that a subject with entirely normal head-motion metrics
carried a large coherent orientation artifact detectable only through
eigenvector-based screening, and argue that motion QC alone is insufficient for
datasets lacking curator-side diffusion quality control.

**Keywords:** diffusion tensor imaging, fractional anisotropy, schizophrenia,
white matter, quality control, confounding

---

## 1. Introduction

### 1.1 The disconnection hypothesis

No single brain region is reliably damaged across Schizophrenia patients, and post-mortem work
has not identified a consistent lesion. This has motivated the *disconnection
hypothesis*: that the disorder arises from disrupted communication **between**
brain regions rather than from dysfunction within any one of them.

White matter carries that communication. It consists of myelinated axons
bundled into tracts connecting distant cortical and subcortical regions. If the
disconnection account is correct, white matter organisation should differ
measurably between patients and controls.

### 1.2 Measuring white matter with diffusion MRI

Diffusion MRI does not image axons directly. It images the thermal motion of
water molecules, which is constrained by local tissue structure. In free fluid,
water diffuses isotropically. Within a coherent axon bundle, diffusion parallel
to the fibres is relatively unimpeded while perpendicular diffusion is hindered
by membranes and myelin sheaths. Diffusion therefore becomes **anisotropic**,
and the degree of anisotropy carries information about tissue organisation.

Measurement follows the Stejskal–Tanner relation:

$$S = S_0 \exp(-b\, \mathbf{g}^{\mathsf{T}} \mathbf{D} \mathbf{g})$$

where $S$ is the diffusion-weighted signal, $S_0$ the signal without diffusion
weighting, $b$ the diffusion weighting factor, $\mathbf{g}$ a unit vector giving
the applied gradient direction, and $\mathbf{D}$ the 3×3 symmetric diffusion
tensor. The quadratic form $\mathbf{g}^{\mathsf{T}}\mathbf{D}\mathbf{g}$ returns
the apparent diffusivity along $\mathbf{g}$.

Taking logarithms linearises the relation in the six unique elements of
$\mathbf{D}$, so that with measurements along many directions the tensor can be
estimated per voxel by least squares. Eigendecomposition of $\mathbf{D}$ yields
eigenvalues $\lambda_1 \ge \lambda_2 \ge \lambda_3$ and the derived scalars:

| Measure | Definition | Interpretation |
|---|---|---|
| FA | $\sqrt{\tfrac{3}{2}}\,\dfrac{\|\mathbf{D}-\mathrm{MD}\cdot\mathbf{I}\|}{\|\mathbf{D}\|}$ | Directional coherence of diffusion, 0–1 |
| MD | $(\lambda_1+\lambda_2+\lambda_3)/3$ | Overall diffusion magnitude |
| AD | $\lambda_1$ | Diffusivity along the principal axis |
| RD | $(\lambda_2+\lambda_3)/2$ | Diffusivity perpendicular to it |

### 1.3 What FA does and does not measure

FA quantifies the *geometric* coherence of water diffusion. Myelination is one
determinant of that coherence, but so are axon density, axon calibre, membrane
integrity, and the degree to which fibres within a voxel run parallel. FA is
also reduced by fibre crossing, which occurs in a substantial proportion of
white matter voxels for reasons that are purely geometric rather than
pathological.

FA is therefore treated throughout this report as a **proxy** for white matter
organisation, not as a measurement of myelin content. Claims are framed
accordingly.

### 1.4 Hypotheses

**Primary (pre-specified).** FA in the corpus callosum is lower in
schizophrenia than in controls. The corpus callosum was designated in advance
as the single primary test because it is the most consistently replicated site
of reduced FA in the schizophrenia literature, and because a study of this size
can support only one confirmatory test.

**Secondary (exploratory).** Group differences in FA across remaining
white matter regions of the JHU ICBM-DTI-81 atlas, and voxelwise across the
tract skeleton via TBSS. All results from these analyses are reported as
exploratory.

**Supporting.** MD, AD and RD are examined in regions showing FA differences,
to characterise whether any effect is driven by parallel or perpendicular
diffusivity.

---

## 2. Methods

### 2.1 Dataset

Data are from the UCLA Consortium for Neuropsychiatric Phenomics LA5c study,
publicly available as OpenNeuro dataset **ds000030**. The release contains 272
participants across four groups: healthy controls (n = 130), schizophrenia
(n = 50), bipolar disorder (n = 49) and ADHD (n = 43). Only the schizophrenia
and control groups were used; bipolar and ADHD participants were excluded by
design at the outset.

Diffusion data were not described in the dataset README, which documents only
anatomical and functional derivatives. Their presence was confirmed by direct
inspection of the raw S3 release. **The absence of a modality from dataset
documentation should not be taken as evidence of its absence from the data.**

Data were retrieved with `aws s3 cp --no-sign-request` from the OpenNeuro S3
bucket. One schizophrenia participant lacked diffusion data, giving 49
potentially eligible patients.

### 2.2 Cohort construction and the protocol confound

An initial cohort was formed by matching each of the 49 eligible schizophrenia
participants to a control of the same sex and similar age, yielding 49 pairs.

Cross-referencing the BIDS JSON sidecars across all participants then revealed
that acquisition was **not uniform**. Data had been collected on two scanners
(serial numbers 35343 and 35426) under three combinations of repetition time
and slice count. Critically, scanner assignment was strongly associated with
diagnosis: 24 of 49 patients versus 8 of 49 controls had been scanned on
scanner 35426 (χ², **p = 0.0001**).

This constitutes an unrecoverable confound. Scanner hardware differs in noise
characteristics, gradient calibration and reconstruction, all of which affect
measured FA. When scanner assignment correlates with diagnosis, the hypotheses
"patients have lower FA" and "scanner 35426 produces lower FA" predict
identical data and cannot be separated by any post-hoc statistical adjustment,
because the two variables carry the same information.

The cohort was therefore rebuilt, matching **exactly on acquisition protocol
and sex** and approximately on age. Feasibility analysis of the control pool by
protocol × sex cell showed that nine male schizophrenia participants scanned on
35426 had no available control counterpart; these were excluded as
unmatchable. Twelve additional control participants were downloaded to complete
the matching. One participant (sub-50029) was excluded for a truncated
acquisition containing 32 rather than 65 volumes.

This produced a matched cohort of **39 pairs (78 participants)**, with
acquisition protocol balanced across groups (p = 1.000) and age comparable
(p = 0.663).

Pair identity was not retained in the cohort file. The matching served to
balance protocol, sex and age at the group level; all analyses are unpaired by
design, with age, sex and motion entered as covariates.

### 2.3 Exclusions

| Stage | Reason | n excluded | Remaining |
|---|---|---|---|
| Eligibility | No diffusion data | 1 SCZ | 49 SCZ |
| Cohort construction | Truncated acquisition (sub-50029) | 1 | — |
| Cohort construction | No protocol-matched control (9 male SCZ, scanner 35426) | 9 | **78** |
| Motion QC | Exceeded max relative displacement threshold | 4 | 74 |
| Protocol | No control on protocol `35343_tr7_sl50` (sub-50005) | 1 | 73 |
| Orientation QC | Coherent orientation artifact (sub-10855) | 1 | **72** |

**Final analysis cohort: 72 participants (36 schizophrenia, 36 control).**

Motion exclusions: sub-10193, sub-11062 (control); sub-50022, sub-50058
(schizophrenia). All four exceeded the pre-specified maximum
volume-to-volume displacement threshold of 3.0 mm.

sub-50005 was the sole remaining participant on acquisition protocol
`35343_tr7_sl50` after its matched control (sub-10193) failed motion QC.
Exclusion follows the rule established during cohort construction, that
participants without a protocol-matched counterpart in the opposite group are
removed.

sub-10855 is discussed in §3.3.

### 2.4 Cohort characteristics

Cross-validation: sub-10321, 5 white matter voxels

voxel (i,j,k)        FSL FA     FSL MD     FSL L1     FSL L3
------------------------------------------------------------
(59,69,23)           0.6131   0.000646   0.001166   0.000301
(45,23,32)           0.6247   0.000609   0.001097   0.000234
(56,9,28)            0.5860   0.000652   0.001116   0.000255
(32,33,23)           0.6395   0.000698   0.001279   0.000264
(44,60,36)           0.5675   0.000771   0.001259   0.000272



### Section 2.4 — Cohort characteristics

| | Control (n = 36) | Schizophrenia (n = 36) | Test |
|---|---|---|---|
| Age, years | 34.167 ± 8.327 | 35.833 ± 9.831 | p = 0.440 |
| Sex, F / M | 10 / 26 | 10 / 26 | p = 1.000 |
| Scanner 35343 / 35426 | 22 / 14 | 22 / 14 | p = 1.000 |
| T1 ghost absent / present | 29 / 7 | 21 / 15 | p = 0.072 |

Scanner, protocol and sex are exactly balanced. The T1 ghosting flag remains
imbalanced; its handling is described in §2.6.


### 2.5 Acquisition parameters

All 78 participants entering preprocessing shared identical relevant
acquisition parameters:

- 65 volumes: 64 diffusion-weighted directions at b = 1000 s/mm² plus one
  b = 0 volume
- Phase encoding direction `j-`, uniform across the cohort
- Total readout time 0.088318 s, uniform across the cohort
- 60 slices (final cohort protocols)

Because phase encoding was uniform and no reverse-encoded acquisition or
fieldmap exists anywhere in the dataset, **susceptibility-induced distortion
cannot be corrected**. This was verified by three independent checks of the
sidecar metadata and directory structure. Consequences are discussed in §6.

A single `acqp.txt` (`0 -1 0 0.088318`) and `index.txt` (65 entries) therefore
served the entire cohort.

### 2.6 Anatomical images excluded from the pipeline

The dataset README documents an aliasing artifact in approximately 20% of T1-
weighted images, appearing as a ghost that may overlap the cortex through one
or both temporal lobes and possibly attributable to a headset worn during
scanning. The prevalence of this artifact is imbalanced across diagnostic
groups: 40% in schizophrenia versus 15% in controls, with every other
diagnostic group at approximately 16%.

Registration algorithms cannot distinguish artifact from anatomy, so any
pipeline routing FA through the participant's T1 would inherit a
group-imbalanced source of misregistration. Registration is therefore performed
**directly from FA to the FMRIB58_FA template**, and T1 images were never
downloaded. This is also the route recommended for TBSS, so the artifact is
avoided at no methodological cost.

### 2.7 Preprocessing

Preprocessing was performed with FSL `[VERIFY version]` and DIPY 1.12.1 on
Apple Silicon hardware. Per participant:

1. **Denoising** — Marchenko–Pastur PCA (DIPY `patch2self` / `mppca`
   `[VERIFY which]`). MP-PCA exploits the redundancy of the 65 volumes: true
   signal occupies few principal components while thermal noise occupies all.
   Marchenko–Pastur theory specifies the eigenvalue distribution expected under
   pure noise, so the signal/noise cut-off is determined analytically rather
   than by a tuning parameter. ~101 s per participant.

2. **Gibbs ringing removal** — Fourier truncation at sharp intensity
   boundaries produces oscillatory artifacts that corrupt FA most severely at
   tissue interfaces. Applied after denoising, since denoising assumes
   independent noise whereas ringing is structured. ~35 s per participant.

3. **Brain extraction** — FSL `bet` applied to the b = 0 volume, producing
   `nodif_brain_mask`.

4. **Eddy current and motion correction** — FSL `eddy_cpu` with `--repol`.
   Eddy jointly models eddy-current-induced distortion (which differs per
   gradient direction) and subject head motion, using a Gaussian process
   prediction of what each volume should look like given the others.
   Slices identified as outliers are replaced by their predictions. ~38 min per
   participant; the cohort was processed at 5-way parallelism.

### 2.8 Tensor fitting

Tensors were fitted with FSL `dtifit` using ordinary least squares, producing
FA, MD, S0, mode, eigenvalues (L1–L3) and eigenvectors (V1–V3). AD was taken as
L1; RD was computed as (L2 + L3)/2 using `fslmaths`, since `dtifit` does not
output these directly.

**Rotated gradient directions were used throughout.** Eddy's motion correction
rotates image volumes, invalidating the original gradient table; eddy therefore
writes `eddy_rotated_bvecs`, and these were supplied to `dtifit`. Use of the
original gradient table produces a normal-appearing FA map with systematically
incorrect orientations, a failure mode that produces no error and is not
detectable from FA alone. The fitting script was written to abort rather than
fall back to the original table if the rotated file was absent.

`[TODO] If the --wls sensitivity analysis is adopted, amend this section and
report both fits.`

### 2.9 Quality control

No curator-side diffusion quality control exists for this dataset. The
published QC protocol (MRIQC) covers anatomical and functional data only, and
the `dwi` flag in `participants.tsv` records scan availability rather than scan
quality. All diffusion QC reported here was performed for this study.

**Stage 1 — Motion.** Eddy's per-volume displacement estimates were extracted
for every participant. Exclusion thresholds were fixed **before any FA value
was inspected**:

- mean volume-to-volume RMS displacement > 2.0 mm
- maximum volume-to-volume RMS displacement > 3.0 mm
- eddy-flagged outlier slices > 10%

**Stage 2 — Orientation.** Colour-encoded FA maps (|V1| modulated by FA) were
generated for every participant and inspected as a montage, then individually
for any participant appearing atypical. The criterion is anatomical: the corpus
callosum must appear red (left–right), the corticospinal tract blue
(superior–inferior), and the superior longitudinal fasciculus green
(anterior–posterior).

Visual inspection was then quantified. For each participant, the proportion of
white matter voxels (FA > 0.25) whose principal eigenvector was strongly
left–right oriented (|V1ₓ| > 0.8) was computed, and an outlier threshold of
median + 3 × IQR applied across the cohort.

### 2.10 Planned analysis

`[TODO — confirm and expand once Phases 6–8 are run]`

**ROI analysis.** Mean FA within each region of the JHU ICBM-DTI-81 atlas,
after nonlinear registration of each participant's FA to FMRIB58_FA space.
Voxels with FA > 1 are excluded from all regional averages (see §3.4).

**Voxelwise analysis.** TBSS. Participant FA maps are registered to
FMRIB58_FA, averaged, and skeletonised to a one-voxel-wide representation of
tract centres, onto which each participant's FA is projected. Skeletonisation
discards voxels near tract edges, where imperfect registration is most damaging
and where FA gradients are steepest.

**Statistics.** Group comparison by general linear model with diagnosis as
the predictor of interest and age, sex and mean relative displacement as
covariates. Voxelwise inference by permutation (FSL `randomise`, 5000
permutations, threshold-free cluster enhancement). Note that `randomise`
outputs 1 − p, so the significance threshold is 0.95.

Effect sizes and 95% confidence intervals are reported alongside p-values
throughout, for reasons given in §6.2.

---

## 3. Quality control results

### 3.1 Preprocessing completion

All 78 participants completed preprocessing. All 74 participants passing motion
QC completed tensor fitting without error.

### 3.2 Head motion

Motion did not differ significantly between groups.


Cross-validation: sub-10321, 5 white matter voxels

voxel (i,j,k)        FSL FA     FSL MD     FSL L1     FSL L3
------------------------------------------------------------
(59,69,23)           0.6131   0.000646   0.001166   0.000301
(45,23,32)           0.6247   0.000609   0.001097   0.000234
(56,9,28)            0.5860   0.000652   0.001116   0.000255
(32,33,23)           0.6395   0.000698   0.001279   0.000264
(44,60,36)           0.5675   0.000771   0.001259   0.000272

| | Control | Schizophrenia |
|---|---|---|
| Mean relative RMS, mm | 0.213 ± 0.075 | 0.238 ± 0.108 |
| &nbsp;&nbsp;median | 0.194 | 0.205 |
| Maximum relative RMS, mm | 0.649 ± 0.345 | 0.816 ± 0.617 |
| &nbsp;&nbsp;median | 0.542 | 0.600 |
| Outlier slices, % | 0.363 ± 0.313 | 0.419 ± 0.397 |
| &nbsp;&nbsp;median | 0.256 | 0.282 |

Welch t(62) = -1.14, p = 0.259; Mann–Whitney p = 0.570; Cohen's d = -0.27, 95% CI [-0.73, 0.20].

The apparent group difference in *maximum* displacement is driven by a small
number of extreme values rather than a group tendency: group medians differ by
0.045 mm while means differ by 0.359 mm, and the schizophrenia standard
deviation is more than double that of controls. A single participant
(sub-50058, maximum 9.04 mm) accounts for much of the difference.

Four participants exceeded the maximum-displacement threshold and were
excluded.

Motion was retained as a GLM covariate despite the absence of a significant
group difference. At this sample size the confidence interval on d is
consistent with a moderate effect, and absence of evidence is not evidence of
absence.



### Section 3.3 — Orientation metric, analysis cohort

- CONTROL 18.875 ± 2.770%, SCZ 19.349 ± 2.570% (p = 0.454, d = -0.18 [-0.64, 0.29])
- Cohort median 19.1%, range 14.0–27.1%
- After excluding sub-10855 the groups are indistinguishable on this metric.


### 3.3 Orientation artifact in sub-10855

Colour-FA inspection identified one participant (sub-10855, control) showing
extensive red signal through frontal and occipital regions, indicating a broad
sheet of apparently coherent left–right diffusion with no anatomical
counterpart.

Quantification confirmed the impression. Across the cohort, the proportion of
white matter voxels strongly left–right oriented had median 19.1% and IQR 3.4%,
giving an outlier threshold of 29.2%. sub-10855 measured **46.4%**, nearly 20
percentage points above the next highest participant (27.1%) and the only
participant exceeding threshold. Mean diffusivity was also elevated
(0.00111 vs a cohort typical value near 0.00093 mm²/s). Eddy's outlier report
showed the first slice flagged repeatedly across volumes at approximately −4
standard deviations, a systematic rather than random pattern.

**This participant's head motion was entirely normal** (mean relative
displacement 0.304 mm, maximum 0.966 mm, 2.4% outlier slices), passing every
motion threshold comfortably.

sub-10855 was excluded. No other participant approached threshold, and the
metric showed no association with diagnosis (control 19.62 ± 5.28,
schizophrenia 19.23 ± 2.64; the elevated control standard deviation is
attributable entirely to this participant).

### 3.4 Degenerate tensor fits

Fractional anisotropy is bounded in [0, 1], yet the maximum FA in every
participant was exactly 1.2247, that is $\sqrt{3/2}$. This value arises when
fitted eigenvalues sum to zero, which occurs when ordinary least squares
returns negative eigenvalues — mathematically permitted by the unconstrained
fit but physically impossible.

Affected voxels comprised 0.27% of in-mask voxels on average (range
0.199–0.364%), with approximately 2% of in-mask voxels having a negative
smallest eigenvalue. None fell in CSF-range mean diffusivity; approximately 36%
fell within white-matter-range mean diffusivity, so these voxels can in
principle fall inside an atlas region.

Voxels with FA > 1 are therefore excluded from all regional averages. This rule
was fixed before any regional FA value was computed.

`[TODO] Report the --wls sensitivity analysis: whether weighted least squares
eliminates the degenerate fits, and whether primary results are unchanged
under it.`

---

## 4. Results

`[TODO — this entire section awaits Phases 6–8]`

### 4.1 Primary hypothesis: corpus callosum

`[TODO] Mean FA in the corpus callosum by group, difference, Cohen's d with
95% CI, p-value from the covariate-adjusted GLM. Report the effect size
prominently; do not lead with the p-value.`

### 4.2 Exploratory regional analysis

`[TODO] Table of JHU atlas regions with group means, effect sizes and
uncorrected and corrected p-values. State the multiple-comparison correction
used. Label the entire section exploratory.`

### 4.3 Exploratory voxelwise analysis

`[TODO] TBSS results. Report skeleton extent, number of significant voxels at
1 − p > 0.95, anatomical location of any clusters. Include a figure of the
skeleton overlaid on the mean FA map.`

### 4.4 Supporting diffusivity measures

`[TODO] MD, AD and RD in any region showing an FA difference, to indicate
whether an effect is driven by parallel or perpendicular diffusivity. Frame
cautiously — see §6.1.`

---

## 5. Discussion

`[TODO — write after Section 4]`

Points to cover:

- Relation of the primary result to the existing schizophrenia DTI literature
  and to the disconnection hypothesis.
- If the primary result is null: what the confidence interval does and does not
  exclude. An underpowered null is weak evidence of absence and must not be
  reported as demonstrating no effect (§6.2).
- Whether any AD/RD pattern is informative, bearing in mind that these measures
  are no more biologically specific than FA.
- The QC finding in §3.3 as a methodological contribution: motion metrics alone
  were insufficient to identify a badly artifacted scan.
- Value of the protocol-matching procedure, and the general point that a
  confound correlated with the grouping variable is invisible in results and
  must be addressed before analysis.

---

## 6. Limitations

### 6.1 Measurement

**No susceptibility distortion correction.** The dataset contains no
reverse-phase-encoded acquisition and no fieldmap, so EPI geometric distortion
could not be corrected. Distortion is most severe near air–tissue interfaces,
affecting orbitofrontal and anterior temporal regions in particular. Because
phase encoding was uniform across all participants, this cannot produce a
spurious group difference, but it renders geometry unreliable in affected
regions and findings there should not be emphasised.

**FA is not a measure of myelin.** FA reflects the directional coherence of
diffusion, which is influenced by axon density, calibre, membrane integrity,
fibre packing and crossing geometry in addition to myelination. Interpretation
in terms of any single microstructural property is not supported by the
measurement.

**Crossing fibres.** The single-tensor model assigns one principal direction
per voxel. In voxels containing crossing fibres — a substantial fraction of
white matter — FA is reduced for geometric reasons unrelated to tissue
integrity.

**No curator-side diffusion QC.** Published quality control for this dataset
covers anatomical and functional data only. All diffusion QC is that described
in §2.9.

### 6.2 Statistical power

With 36 participants per group, power to detect an effect of the magnitude
typically reported in this literature (d ≈ 0.4) is approximately 40%. Detecting
an effect with 80% power would require d ≈ 0.64.

Three consequences follow, and are applied throughout:

1. A single primary hypothesis was pre-specified; all other analyses are
   exploratory.
2. Effect sizes and confidence intervals are reported alongside p-values.
3. A null result constitutes weak evidence of absence only. It constrains the
   plausible effect size and contributes to future meta-analysis, but cannot
   establish that no difference exists.

### 6.3 Sample

Participants excluded for protocol mismatch were disproportionately male
patients scanned on one scanner, so the final sample is not fully
representative of the original patient group. Clinical variables including
medication status, illness duration and symptom severity were not modelled.

The T1 ghosting flag remains imbalanced across groups (p = 0.072). Because
anatomical images were excluded from the pipeline entirely, this does not
affect the present analysis, but it may indicate systematic differences in how
the two groups were scanned or positioned.

---

## 7. Reproducibility

All processing scripts, the cohort definition file, and a complete log of
methodological decisions with their evidentiary basis are retained in the
project repository.

| File | Contents |
|---|---|
| `cohort_final.tsv` | Cohort definition with acquisition metadata |
| `cohort_decisions.md` | Every methodological decision, its evidence and consequences |
| `derivatives/motion_qc.tsv` | Per-participant motion metrics and exclusion flags |
| `derivatives/orientation_qc.tsv` | Per-participant orientation metrics |
| `code/process_one_subject.sh` | Per-participant preprocessing |
| `code/run_batch.sh` | Parallel batch driver |
| `code/phase4f_motion_qc.py` | Motion QC and exclusion |
| `code/phase5_dtifit.sh` | Tensor fitting |
| `code/phase5_colorfa_check.py` | Orientation QC |

All analysis subject lists are derived from `derivatives/motion_qc.tsv`
filtered on `exclude_final`, never from directory listings, since the
`raw/` and `derivatives/dti/` trees contain participants not in the analysis
cohort.

---

## References

`[TODO — complete]`

1. Poldrack RA, et al. A phenome-wide examination of neural and cognitive
   function. *Scientific Data* 2016. `[dataset reference — verify]`
2. Basser PJ, Mattiello J, LeBihan D. MR diffusion tensor spectroscopy and
   imaging. *Biophysical Journal* 1994.
3. Stejskal EO, Tanner JE. Spin diffusion measurements. *Journal of Chemical
   Physics* 1965.
4. Smith SM, et al. Tract-based spatial statistics. *NeuroImage* 2006.
5. Andersson JLR, Sotiropoulos SN. An integrated approach to correction for
   off-resonance effects and subject movement in diffusion MR imaging.
   *NeuroImage* 2016.
6. Veraart J, et al. Denoising of diffusion MRI using random matrix theory.
   *NeuroImage* 2016.
7. Kellner E, et al. Gibbs-ringing artifact removal based on local subvoxel
   shifts. *Magnetic Resonance in Medicine* 2016.
8. Mori S, et al. JHU ICBM-DTI-81 white matter atlas. `[verify citation]`
9. Winkler AM, et al. Permutation inference for the general linear model.
   *NeuroImage* 2014.
10. `[TODO]` A schizophrenia DTI meta-analysis supporting the corpus callosum
    primary hypothesis — e.g. Kelly S, et al. ENIGMA Schizophrenia DTI Working
    Group. *Molecular Psychiatry* 2018. `[verify]`

---

## Appendix A — Decision log summary

| Decision | Evidence | Consequence |
|---|---|---|
| Restrict to SCZ vs control | Study scope | Bipolar and ADHD groups unused |
| Match on acquisition protocol | Scanner × diagnosis χ² p = 0.0001 | 9 patients excluded as unmatchable |
| Exclude sub-50029 | 32 volumes rather than 65 | — |
| Exclude T1 from pipeline | Ghosting 40% SCZ vs 15% control | FA registered directly to FMRIB58_FA |
| No susceptibility correction | No reverse-PE run, no fieldmap | Frontal/temporal geometry unreliable |
| Motion thresholds fixed before FA inspection | — | 4 exclusions |
| Exclude sub-50005 | Sole participant on `tr7_sl50` after partner failed QC | Consistent with unmatchable rule |
| Exclude sub-10855 | L–R orientation 46.4% vs median 19.1%, threshold 29.2% | Motion was normal; only orientation QC caught it |
| Exclude FA > 1 voxels from ROI means | 0.27% of voxels, $\sqrt{3/2}$ degenerate-fit signature | Applied to all regional averages |
| Unpaired analysis | Pair identity not retained | Age, sex, motion as covariates |
| Corpus callosum pre-specified as primary | ~40% power; one confirmatory test only | All else exploratory |