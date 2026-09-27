# White Matter Microstructure in Schizophrenia: A Diffusion Tensor Imaging Study of the UCLA Consortium for Neuropsychiatric Phenomics Cohort

**Frederick Lancia**
Erasmus Mundus Maths DISC — University of Verona
Diffusion MRI Mini-Course Project

---

## Abstract

Reduced fractional anisotropy (FA) in white matter is among the most replicated
neuroimaging findings in schizophrenia, and is commonly interpreted as evidence
for the disconnection hypothesis [1], which holds that psychiatric and
neurological disorders stem from a breakdown in the functional integration and
communication between distributed brain regions rather than damage to isolated
centres. We tested for group differences in white matter FA between individuals
with schizophrenia and healthy controls using diffusion-weighted imaging from
the UCLA Consortium for Neuropsychiatric Phenomics LA5c dataset (OpenNeuro
ds000030) [2].

From 272 available participants we constructed a cohort matched exactly on
acquisition protocol and sex and approximately on age. This matching was
necessary because acquisition protocol was strongly confounded with diagnosis
in the full sample (χ², p = 0.0001), an imbalance that would have rendered any
FA finding uninterpretable. After preprocessing and a two-stage quality control
procedure, 72 participants (36 schizophrenia, 36 control) entered analysis.

The pre-specified primary hypothesis was not supported: corpus callosum FA did
not differ significantly between groups (adjusted difference −0.0064, 95% CI
[−0.0162, +0.0034], d = −0.29, p = 0.197). Exploratory analysis instead
indicated a **whole-brain** reduction in FA in schizophrenia (mean FA across the
JHU atlas d = −0.58, p = 0.018), with lower FA in 45 of 50 regions and nine
regions surviving FDR correction. Voxelwise analysis (TBSS, 5000 permutations,
TFCE-corrected) corroborated this, with 13.9% of the white matter skeleton
showing lower FA in schizophrenia, distributed across 34 of 49 atlas tracts,
and no voxel showing the reverse. The pattern was driven predominantly by
increased radial rather than decreased axial diffusivity, and was not
attributable to head motion, scan-quality differences, or the dataset's
documented imaging artifact. These results closely parallel the largest
coordinated analysis in this literature [3].

We additionally report that a participant with entirely normal head-motion
metrics carried a large coherent orientation artifact detectable only through
eigenvector-based screening, and argue that motion QC alone is insufficient for
datasets lacking curator-side diffusion quality control.

**Keywords:** diffusion tensor imaging, fractional anisotropy, schizophrenia,
white matter, quality control, confounding

---

## 1. Introduction

### 1.1 The disconnection hypothesis

No single brain region is reliably damaged across schizophrenia patients, and
post-mortem work has not identified a consistent lesion. This has motivated the
*disconnection hypothesis* [1, 4]: that the disorder arises from disrupted
communication **between** brain regions rather than from dysfunction within any
one of them.

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

Measurement follows the Stejskal–Tanner relation [5]:

$$S = S_0 \exp(-b\, \mathbf{g}^{\mathsf{T}} \mathbf{D} \mathbf{g})$$

where $S$ is the diffusion-weighted signal, $S_0$ the signal without diffusion
weighting, $b$ the diffusion weighting factor, $\mathbf{g}$ a unit vector giving
the applied gradient direction, and $\mathbf{D}$ the 3×3 symmetric diffusion
tensor. The quadratic form $\mathbf{g}^{\mathsf{T}}\mathbf{D}\mathbf{g}$ returns
the apparent diffusivity along $\mathbf{g}$.

Taking logarithms linearises the relation in the six unique elements of
$\mathbf{D}$, so that with measurements along many directions the tensor can be
estimated per voxel by least squares [6]. Eigendecomposition of $\mathbf{D}$
yields eigenvalues $\lambda_1 \ge \lambda_2 \ge \lambda_3$ and the derived
scalars:

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
as the single primary test because it is among the most consistently replicated
sites of reduced FA in the schizophrenia literature [3], and because a study of
this size can support only one confirmatory test.

**Secondary (exploratory).** Group differences in FA across remaining
white matter regions of the JHU ICBM-DTI-81 atlas [7–9], and voxelwise across
the tract skeleton via TBSS [10]. All results from these analyses are reported
as exploratory.

**Supporting.** MD, AD and RD are examined in regions showing FA differences,
to characterise whether any effect is driven by parallel or perpendicular
diffusivity.

---

## 2. Methods

### 2.1 Dataset

Data are from the UCLA Consortium for Neuropsychiatric Phenomics LA5c study,
publicly available as OpenNeuro dataset **ds000030** [2]. The release contains
272 participants across four groups: healthy controls (n = 130), schizophrenia
(n = 50), bipolar disorder (n = 49) and ADHD (n = 43). Only the schizophrenia
and control groups were used; bipolar and ADHD participants were excluded by
design at the outset.

Diffusion data were not described in the dataset README, which documents only
anatomical and functional derivatives. Their presence was confirmed by direct
inspection of the raw S3 release.

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

| | Control (n = 36) | Schizophrenia (n = 36) | Test |
|---|---|---|---|
| Age, years | 34.17 ± 8.33 | 35.83 ± 9.83 | p = 0.440 |
| Sex, F / M | 10 / 26 | 10 / 26 | p = 1.000 |
| Scanner 35343 / 35426 | 22 / 14 | 22 / 14 | p = 1.000 |
| Protocol `tr8.4_sl60` / `tr9_sl60` | 22 / 14 | 22 / 14 | p = 1.000 |
| T1 ghost absent / present | 29 / 7 | 21 / 15 | p = 0.072 |

Scanner, protocol and sex are exactly balanced. The T1 ghosting flag remains
imbalanced; its handling is described in §2.6 and its residual implications in
§6.3.

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
sidecar metadata and directory structure. Consequences are discussed in §6.1.

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
downloaded. This is also the route recommended for TBSS [10], so the artifact
is avoided at no methodological cost.

### 2.7 Preprocessing

Preprocessing used FSL 6.0.7.21 [11] and DIPY 1.12.1 [12] on Apple
Silicon hardware, with nibabel 5.4.2, numpy 2.5.2, scipy 1.18.0, pandas 3.0.5
and statsmodels 0.14.6. Per participant:

1. **Denoising** — Marchenko–Pastur PCA (DIPY `mppca`) [13]. MP-PCA exploits
   the redundancy of the 65 volumes: true signal occupies few principal
   components while thermal noise occupies all. Marchenko–Pastur theory
   specifies the eigenvalue distribution expected under pure noise, so the
   signal/noise cut-off is determined analytically rather than by a tuning
   parameter. ~101 s per participant.

2. **Gibbs ringing removal** — DIPY `gibbs_removal` [14]. Fourier truncation at
   sharp intensity boundaries produces oscillatory artifacts that corrupt FA
   most severely at tissue interfaces. Applied after denoising, since denoising
   assumes independent noise whereas ringing is structured. ~35 s per
   participant.

3. **Brain extraction** — FSL `bet` [15] applied to the b = 0 volume, producing
   `nodif_brain_mask`.

4. **Eddy current and motion correction** — FSL `eddy_cpu` [16] with `--repol`
   [17]. Eddy jointly models eddy-current-induced distortion (which differs per
   gradient direction) and subject head motion, using a Gaussian process
   prediction of what each volume should look like given the others. Slices
   identified as outliers are replaced by their predictions. ~38 min per
   participant; the cohort was processed at 5-way parallelism.

### 2.8 Tensor fitting

Tensors were fitted with FSL `dtifit` using ordinary least squares, producing
FA, MD, S0, mode, eigenvalues (L1–L3) and eigenvectors (V1–V3). AD was taken as
L1; RD was computed as (L2 + L3)/2 using `fslmaths`, since `dtifit` does not
output these directly.

**Rotated gradient directions were used throughout.** A diffusion measurement
is a pair: an image volume and the gradient direction along which it was
sensitised, the latter expressed in the scanner's coordinate frame. When eddy
resamples a volume to correct head rotation, the correspondence between the
voxel grid and that direction changes, so the gradient table must be rotated by
the same amount. Eddy writes `eddy_rotated_bvecs` for this purpose, and these
were supplied to `dtifit`.

Using the original gradient table produces a normal-appearing FA map. Because
rotation preserves eigenvalues, a uniform rotation error leaves FA, MD, AD and
RD numerically unchanged and corrupts only the eigenvector directions. Because
head rotation differs between volumes, the errors are mutually inconsistent
across the 65 measurements, and the least-squares fit resolves that
inconsistency by returning a tensor biased toward isotropy — depressing FA in
proportion to how much the participant moved. As patients move more on average,
this failure mode can manufacture a group difference from a filename error,
while producing no error message and no anomaly in any scalar map. The fitting
script was therefore written to abort rather than fall back to the original
table if the rotated file was absent, and orientation was verified explicitly
(§2.9, §3.3).

A weighted least squares alternative (`dtifit --wls`) was evaluated and
rejected; see §3.4.

### 2.9 Quality control

No curator-side diffusion quality control exists for this dataset. The
published QC protocol (MRIQC) covers anatomical and functional data only, and
the `dwi` flag in `participants.tsv` records scan availability rather than scan
quality. All diffusion QC reported here was performed for this study.

**Stage 1 — Motion.** Eddy's per-volume displacement estimates were extracted
for every participant. A participant was excluded if **any** of the following
held:

- mean volume-to-volume RMS displacement exceeded 2.0 mm
- maximum volume-to-volume RMS displacement exceeded 3.0 mm
- more than 10% of slices were flagged as outliers by eddy

These thresholds were fixed **before any FA value was inspected**. Setting a
threshold after observing which participants have low FA would allow the
criterion to be chosen so as to remove inconvenient data.

**Stage 2 — Orientation.** Colour-encoded FA maps (|V1| modulated by FA) were
generated for every participant and inspected as a montage, then individually
for any participant appearing atypical. The criterion is anatomical: the corpus
callosum must appear red (left–right), the corticospinal tract blue
(superior–inferior), and the superior longitudinal fasciculus green
(anterior–posterior).

Visual inspection was then quantified. For each participant, the proportion of
white matter voxels (FA > 0.25) whose principal eigenvector was strongly
left–right oriented (|V1ₓ| > 0.8) was computed, and a participant was excluded
if this proportion exceeded the cohort median plus three interquartile ranges.

### 2.10 Registration and analysis

**Registration.** Every brain differs in size, shape and position in the
scanner, so a given voxel index corresponds to different anatomy in different
participants; group comparison requires that the same index mean the same
anatomy. Each participant's FA map was therefore warped onto the FMRIB58_FA
1 mm template, an average FA map derived from 58 individuals. This proceeded in
two stages: a 12-degree-of-freedom affine alignment (FSL `flirt`) [18], which
fixes global position, scale and shear but cannot accommodate local shape
differences, followed by nonlinear warping (FSL `fnirt`, configuration
`FA_2_FMRIB58_1mm`) [19], which assigns a separate displacement to each voxel
and so aligns local structure. The resulting warp was applied to that
participant's MD, AD and RD maps. No anatomical image was used at any stage.
Registration quality is reported in §3.5.

**ROI analysis.** Mean FA, MD, AD and RD within each region of the JHU
ICBM-DTI-81 atlas [7–9], a set of hand-segmented white matter labels defined in
template space. Voxels with FA > 1 were excluded from all regional averages
(§3.4), as were regions containing fewer than 10 valid voxels.

**Voxelwise analysis.** TBSS [10]. Participant FA maps are registered to
FMRIB58_FA, averaged, and skeletonised at FA > 0.2 to a one-voxel-wide
representation of tract centres, onto which each participant's FA is projected.
Skeletonisation discards voxels near tract edges, where imperfect registration
is most damaging and where FA gradients are steepest.

**Statistics.** Regional group comparison by ordinary least squares regression
of the form `measure ~ diagnosis + age + sex + mean relative displacement`,
corrected across the 50 atlas regions by the Benjamini–Hochberg false discovery
rate at q < 0.05 [20]. Voxelwise inference by permutation (FSL `randomise`,
5000 permutations) [21] with threshold-free cluster enhancement [22], using the
same covariates, demeaned, and two contrasts (control > schizophrenia and its
reverse); note that `randomise` outputs 1 − p, so the significance threshold is
0.95. Because `randomise` matches design-matrix rows to participants by
position, the design matrix and the TBSS input order were generated from a
single sorted list and verified to match before inference was run.

Effect sizes and 95% confidence intervals are reported alongside p-values
throughout, for reasons given in §6.2. Reported Cohen's d values are computed
from unadjusted group means, while p-values derive from the covariate-adjusted
model.

---

## 3. Quality control results

### 3.1 Processing completion

No participant was lost to processing failure at any stage. All 78 participants
entering preprocessing completed it; all 74 remaining after motion exclusions
completed tensor fitting; all 72 remaining after the protocol and orientation
exclusions completed registration. The cohort reduced from 78 to 72 solely
through the exclusion criteria in §2.3.

### 3.2 Head motion

Motion did not differ significantly between groups.

| | Control | Schizophrenia |
|---|---|---|
| Mean relative RMS, mm | 0.213 ± 0.075 | 0.238 ± 0.108 |
| &nbsp;&nbsp;median | 0.194 | 0.205 |
| Maximum relative RMS, mm | 0.649 ± 0.345 | 0.816 ± 0.617 |
| &nbsp;&nbsp;median | 0.542 | 0.600 |
| Outlier slices, % | 0.363 ± 0.313 | 0.419 ± 0.397 |
| &nbsp;&nbsp;median | 0.256 | 0.282 |

Welch t(62) = −1.14, p = 0.259; Mann–Whitney p = 0.570; Cohen's d = −0.27,
95% CI [−0.73, 0.20].

At the pre-exclusion stage (n = 78) the apparent group difference in *maximum*
displacement was driven by a small number of extreme values rather than a group
tendency: group medians differed by 0.045 mm while means differed by 0.359 mm,
and the schizophrenia standard deviation was more than double that of controls.
A single participant (sub-50058, maximum 9.04 mm) accounted for much of the
difference. Four participants exceeded the maximum-displacement threshold and
were excluded.

Motion was retained as a GLM covariate despite the absence of a significant
group difference, since the confidence interval on d remains consistent with a
moderate effect.

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

sub-10855 was excluded. Within the resulting analysis cohort the metric showed
no association with diagnosis (control 18.88 ± 2.77%, schizophrenia
19.35 ± 2.57%; p = 0.454, d = −0.18, 95% CI [−0.64, 0.29]; cohort median 19.1%,
range 14.0–27.1%).

Four further participants flagged on visual inspection (sub-10460, sub-10530,
sub-50007, sub-50032) were examined at three slice levels and found normal.
Their atypical appearance in the montage arose because the montage samples a
single slice at 50% of each participant's slice stack; since participants
differ in slice count and head position, that slice falls at different
anatomical levels in different people, and a brain sampled low appears smaller
and differently shaped. A montage is adequate for detecting gross orientation
failure and inadequate for judging anything subtler.

### 3.4 Degenerate tensor fits

Fractional anisotropy is bounded in [0, 1], yet the maximum FA in every
participant was exactly 1.2247, that is $\sqrt{3/2}$.

The explanation is that the fit is not constrained to be physically possible.
Ordinary least squares returns whichever six tensor elements minimise squared
error, with no requirement that the resulting eigenvalues be positive. In
voxels where the signal approaches the noise floor, the measurements are
mutually inconsistent and the best-fitting tensor can include negative
eigenvalues, implying negative diffusivity along some axis. Because FA divides
the spread of the eigenvalues by their overall magnitude, negative eigenvalues
cancel positive ones in the denominator and inflate the ratio; $\sqrt{3/2}$ is
the mathematical ceiling, attained when the eigenvalues sum to zero. Identical
maxima across participants therefore reflect a property of the estimator, not a
measurement.

Approximately 2% of in-mask voxels (voxels within the BET brain mask) had a
negative smallest eigenvalue — the broad symptom. The severe subset, FA > 1,
comprised 0.27% of in-mask voxels on average (range 0.199–0.364%).

These voxels were characterised to determine whether they could affect the
analysis. Mean diffusivity distinguishes tissue types: CSF approaches free
water at roughly 0.003 mm²/s, white matter roughly 0.0007–0.0009 mm²/s. **None**
of the affected voxels fell in the CSF range, so they are not simply free-water
fits in the ventricles. Approximately 36% fell in the white-matter range, and
approximately 40% survived eroding the brain mask by three voxels, so they are
neither confined to the brain periphery nor to non-tissue, and can fall inside
an atlas region.

A weighted least squares sensitivity analysis (`dtifit --wls`, three
participants) **increased** rather than reduced the number of degenerate fits
(OLS 710 / 671 / 575 voxels versus WLS 746 / 710 / 833). FSL derives WLS
weights from an initial OLS fit, so in noise-dominated voxels the weights are
themselves corrupted and amplify the affected measurements. Ordinary least
squares was therefore retained.

Voxels with FA > 1 are excluded from all regional averages. This rule was fixed
before any regional FA value was computed.

### 3.5 Pipeline validation

**Independent reimplementation.** The tensor fit was independently implemented
from first principles (design matrix construction from the rotated gradient
table, solution by pseudoinverse, eigendecomposition) and compared against FSL
at five white matter voxels in one participant. The two implementations agreed
to four decimal places on both FA and MD:

| Voxel (i,j,k) | FSL FA | Independent FA | FSL MD | Independent MD |
|---|---|---|---|---|
| (59,69,23) | 0.6131 | 0.6131 | 0.000646 | 0.000646 |
| (45,23,32) | 0.6247 | 0.6247 | 0.000609 | 0.000609 |
| (56,9,28) | 0.5860 | 0.5860 | 0.000652 | 0.000652 |
| (32,33,23) | 0.6395 | 0.6395 | 0.000698 | 0.000698 |
| (44,60,36) | 0.5675 | 0.5675 | 0.000771 | 0.000771 |

This confirms consistent handling of the gradient table, b-values, log
transform and eigenvalue ordering between the two implementations.

**Registration quality.** All 72 standard-space FA maps were averaged. The
cohort mean FA correlated with the FMRIB58_FA template at **r = 0.974** within
template white matter, with mean FA of 0.346, and showed sharply delineated
tract structure including the corpus callosum, internal capsule and corona
radiata. Misregistration would produce a blurred average lacking recognisable
tract anatomy.

---

## 4. Results

### 4.1 Primary hypothesis: corpus callosum

**The pre-specified primary hypothesis was not supported.**

| | Control | Schizophrenia |
|---|---|---|
| Corpus callosum mean FA | 0.6412 | 0.6352 |

Adjusted difference (diagnosis, controlling for age, sex and motion):
**−0.0064, 95% CI [−0.0162, +0.0034]**, Cohen's d = −0.29, p = 0.197.

The direction of effect is consistent with the hypothesis and the effect size
is small to moderate, but the confidence interval includes zero. The interval
constrains any true difference to approximately 0.016 FA units or less, roughly
2.5% of the regional mean.

Subregional effects were heterogeneous:

| Subregion | Cohen's d | p |
|---|---|---|
| Genu of corpus callosum | −0.44 | 0.059 |
| Body of corpus callosum | −0.17 | 0.437 |
| Splenium of corpus callosum | −0.24 | 0.258 |

The genu shows the largest effect and the body the smallest. Averaging the
three subregions into a single primary ROI therefore diluted the strongest
signal. This is a limitation of the pre-specified ROI definition and is
discussed in §5.2; the primary hypothesis **cannot** be redefined post hoc.

### 4.2 Exploratory regional analysis: a whole-brain FA reduction

*All results in this section are exploratory and hypothesis-generating.*

**The effect is global rather than regionally selective.** Mean FA averaged
across all 50 JHU atlas regions was lower in schizophrenia (control 0.5161,
schizophrenia 0.5070; adjusted difference −0.0090, **d = −0.58, 95% CI [−1.06,
−0.11], p = 0.018**). Lower FA in schizophrenia was observed in **45 of 50
regions**. A sign test on that count gives p = 4.2 × 10⁻⁹, but this figure is
anticonservative: the sign test assumes independent observations, whereas a
participant with slightly lower FA throughout contributes 50 correlated votes
rather than 50 independent ones. The count is reported descriptively, and the
global test above — which collapses the correlated regions into one value per
participant before testing — is the appropriate formal inference.

Nine regions survived FDR correction [20] at q < 0.05:

| Region | Control FA | SCZ FA | d | p | p(FDR) |
|---|---|---|---|---|---|
| Posterior thalamic radiation R | 0.5446 | 0.5217 | −0.91 | 0.0004 | 0.015 |
| Retrolenticular internal capsule L | 0.4393 | 0.4162 | −0.87 | 0.0006 | 0.015 |
| Posterior thalamic radiation L | 0.5224 | 0.5009 | −0.81 | 0.0016 | 0.026 |
| Cingulum (cingulate gyrus) L | 0.4057 | 0.3852 | −0.65 | 0.0026 | 0.026 |
| Sagittal stratum R | 0.4836 | 0.4666 | −0.80 | 0.0026 | 0.026 |
| Anterior corona radiata R | 0.4430 | 0.4230 | −0.74 | 0.0036 | 0.030 |
| Posterior limb of internal capsule L | 0.5705 | 0.5515 | −0.75 | 0.0042 | 0.030 |
| Posterior corona radiata L | 0.5647 | 0.5458 | −0.65 | 0.0080 | 0.049 |
| Inferior fronto-occipital fasciculus R | 0.4375 | 0.4217 | −0.62 | 0.0087 | 0.049 |

Given the global shift, these regions should not be read as selectively
affected; the reasoning is set out in §5.3.

**The effect is not attributable to scan quality.** No quality metric differed
between groups, and none correlated meaningfully with global FA:

| Metric | Group difference | Correlation with global FA |
|---|---|---|
| Mean relative displacement | p = 0.259 | +0.063 |
| Outlier slices, % | p = 0.506 | +0.003 |
| Maximum displacement | p = 0.162 | −0.041 |

A quality confound would require both a group difference in the metric and an
association with FA. Neither is present.

**The effect is not driven by individual participants.** Leave-one-out analysis
of the strongest region (posterior thalamic radiation R, full-sample d = −0.91)
showed a maximum shift of 0.111 in Cohen's d across all 72 exclusions, below
the 0.15 threshold set in advance. Every influential exclusion moved the effect
*away* from zero, so no participant is inflating the result.

Figure: `derivatives/figures_forest_FA.png` (regional effect sizes with 95%
confidence intervals, FDR-significant regions highlighted).

### 4.3 Exploratory voxelwise analysis (TBSS)

*Exploratory.*

The study-specific mean FA skeleton comprised **114,263 voxels**.

| Contrast | Significant voxels (1 − p > 0.95) | % of skeleton | Max 1 − p |
|---|---|---|---|
| Control > schizophrenia | **15,887** | 13.9% | 0.9904 |
| Schizophrenia > control | **0** | 0.0% | 0.4106 |

The result is strongly **directional**. Almost 14% of the skeleton showed lower
FA in schizophrenia after TFCE correction across all 114,263 voxels, while not
a single voxel showed the reverse, and the reverse contrast did not approach
threshold anywhere (peak 1 − p = 0.41, i.e. p = 0.59). Noise and residual
artifact typically produce scattered significance in both directions; a
one-sided result of this extent does not have that character.

The result is also **diffuse rather than focal**. Significant voxels were
distributed across **34 of 49** labelled atlas tracts, with 6,449 voxels (41%
of all significant voxels) falling outside labelled JHU regions entirely. The
most affected tracts, as a proportion of each tract's skeleton:

| Tract | Significant / skeleton voxels | % |
|---|---|---|
| Retrolenticular part of internal capsule L | 1165 / 1599 | 72.9 |
| Inferior fronto-occipital fasciculus R | 301 / 457 | 65.9 |
| Posterior thalamic radiation L | 365 / 583 | 62.6 |
| Anterior corona radiata R | 1094 / 1768 | 61.9 |
| Superior fronto-occipital fasciculus L | 40 / 66 | 60.6 |
| Posterior thalamic radiation R | 623 / 1111 | 56.1 |
| Sagittal stratum R | 239 / 478 | 50.0 |
| Posterior limb of internal capsule L | 370 / 742 | 49.9 |
| Retrolenticular part of internal capsule R | 374 / 767 | 48.8 |
| Pontine crossing tract | 860 / 1768 | 48.6 |
| Cingulum (hippocampus) L | 140 / 289 | 48.4 |
| Cerebral peduncle R | 301 / 632 | 47.6 |
| Posterior corona radiata L | 510 / 1147 | 44.5 |
| Anterior limb of internal capsule R | 345 / 807 | 42.8 |

The voxelwise result **converges with the regional analysis** (§4.2): the
retrolenticular internal capsule, posterior thalamic radiation bilaterally,
anterior corona radiata, sagittal stratum, posterior limb of the internal
capsule, posterior corona radiata and inferior fronto-occipital fasciculus
appear prominently in both, despite the two analyses differing in spatial
resolution, statistical model and correction procedure.

Two observations temper the picture. First, peak significance is modest: the
strongest voxel reaches p = 0.0096, and no voxel approaches p < 0.001. The
finding is broad rather than deep, which is what a moderate global shift
predicts. Second, brainstem structures appear among the most affected tracts
(pontine crossing tract 48.6%, cerebral peduncle R 47.6%), and these lie where
uncorrected susceptibility distortion is most severe (§6.1).

Figure: `derivatives/figures_tbss_tstat1.png` (skeleton in green, significant
voxels in red/yellow, over the study mean FA).

### 4.4 Supporting diffusivity measures

Across the nine FDR-significant regions, the pattern was dominated by
**increased radial diffusivity** rather than decreased axial diffusivity:

| Region | FA d | AD d | RD d | Pattern |
|---|---|---|---|---|
| Posterior thalamic radiation R | −0.91 | −0.39 | +0.23 | both |
| Retrolenticular internal capsule L | −0.87 | −0.19 | +0.57 | RD only |
| Posterior thalamic radiation L | −0.81 | −0.41 | +0.46 | both |
| Cingulum (cingulate gyrus) L | −0.65 | −0.18 | +0.50 | RD only |
| Sagittal stratum R | −0.80 | −0.31 | +0.42 | both |
| Anterior corona radiata R | −0.74 | +0.09 | +0.60 | RD only |
| Posterior limb of internal capsule L | −0.75 | −0.33 | +0.56 | both |
| Posterior corona radiata L | −0.65 | −0.32 | +0.22 | both |
| Inferior fronto-occipital fasciculus R | −0.62 | −0.69 | +0.30 | both |

RD was elevated in all nine regions (d = +0.22 to +0.60), while AD was reduced
in eight of nine but generally more weakly. In three regions the pattern was
RD-only. This matches the direction reported by the ENIGMA Schizophrenia DTI
Working Group, which found significantly higher mean and radial diffusivity in
patients and larger effects for FA than for diffusivity measures [3].

In the animal literature an RD-dominant pattern has been associated with
myelin-related change, but that mapping does not transfer reliably to human
group comparisons at this resolution, where crossing fibres, partial volume and
registration error all influence RD. The pattern is reported descriptively; no
microstructural mechanism is claimed (§6.1).

---

## 5. Discussion

### 5.1 Summary

The pre-specified test of reduced corpus callosum FA in schizophrenia was not
significant (d = −0.29, 95% CI [−0.016, +0.003] in FA units, p = 0.197).
Exploratory analysis instead indicated a **brain-wide** reduction in FA
(d = −0.58, p = 0.018), present in 45 of 50 atlas regions and surviving FDR
correction in nine. Voxelwise analysis corroborated this independently: 13.9%
of the white matter skeleton showed lower FA in schizophrenia after permutation
testing with TFCE correction, spread across 34 of 49 tracts, with no voxel
showing the opposite direction. The effect was not explained by head motion,
scan quality, the dataset's documented imaging artifact, or individual
participants, and was dominated by increased radial diffusivity.

### 5.2 The primary null

Three readings are available, and the data do not cleanly separate them.

**Insufficient power.** With 36 per group, power to detect d = 0.4 is
approximately 40%, and the observed corpus callosum effect (d = −0.29) is
smaller than that. ENIGMA, pooling 4,322 individuals, reported d = 0.39 for the
corpus callosum (body 0.39, genu 0.37) [3] — magnitudes this study was not
powered to detect. A null was therefore the more likely outcome even if the
published effect is real, and this is the most parsimonious reading.

**ROI definition.** The primary ROI pooled genu, body and splenium. The genu
alone reached d = −0.44 (p = 0.059) while the body reached only d = −0.17.
Where an effect occupies part of a region, averaging over the whole region
mixes affected with unaffected tissue and shrinks the measured difference, so a
narrower pre-specification might have been more sensitive. Redefining the
primary ROI after seeing these numbers would convert a confirmatory test into
an exploratory one.

**Genuine relative sparing.** The corpus callosum effect (d = −0.29) is
approximately half the global effect (d = −0.58), and the corpus callosum does
not appear among the most affected tracts in the voxelwise analysis. This would
be consistent with it being *less* affected than white matter generally in this
sample — though it sits against ENIGMA, where the corpus callosum was among the
strongest regions [3]. Given the confidence intervals involved this remains
speculative.

What the result does **not** support is a claim that corpus callosum FA is
unaffected in schizophrenia. The confidence interval excludes large effects and
is entirely consistent with the published small-to-moderate effects.

### 5.3 The global finding

A diffuse rather than focal difference is consistent with the disconnection
hypothesis in its broad form [1, 4]. It also closely parallels ENIGMA [3],
which found significant FA reductions in 20 of 25 regions of interest, a
whole-skeleton effect of d = 0.42, and the largest regional effects in the
anterior corona radiata (d = 0.40) and corpus callosum (d = 0.39). The global
effect observed here (d = −0.58) is somewhat larger but within the range
expected given this sample size, and the anterior corona radiata is prominent
in both analyses.

Four features of the data support taking the finding seriously. The regional
and voxelwise analyses agree despite using different spatial resolutions,
models and correction procedures. The voxelwise result is strictly
one-directional, which residual (post-correction) artifact would not typically
produce. No measured quality metric differed between groups or predicted global
FA. And the effect was robust to removing any individual participant.

Three cautions apply.

First, the finding is **exploratory**. It was not pre-specified, and it emerged
after the primary test returned null. It requires independent replication
before being treated as established.

Second, a global shift makes regional interpretation hazardous. If every tract
is reduced by a similar amount, the tracts that cross a significance threshold
are those where FA can be measured most precisely — large tracts, tracts that
register consistently, tracts with coherent single-fibre architecture. The
FDR-significant regions and the most affected TBSS tracts are therefore better
described as the places where a brain-wide difference was most detectable than
as selectively affected pathways. The marked left–right asymmetry in the
anterior corona radiata (61.9% of the right skeleton significant versus 24.6%
of the left) illustrates this: no plausible mechanism makes a tract 2.5 times
more affected on one side.

Third, a global FA difference is what a residual systematic difference in data
quality would produce. The quality checks (§4.2) and the ghosting sensitivity
analysis (§6.3) argue against this, but cannot exclude a source not captured by
the available metrics.

### 5.4 Methodological contributions

**Protocol confounding must be checked before analysis.** Acquisition protocol
was strongly associated with diagnosis in the full sample (p = 0.0001). Had
this gone unexamined, "patients have lower FA" and "scanner 35426 produces
lower FA" would have predicted identical data, with no post-hoc adjustment
capable of separating them. Detecting it required inspecting BIDS sidecar
metadata rather than the dataset documentation, which described neither the
protocol heterogeneity nor the diffusion data itself.

**Motion QC alone is insufficient.** sub-10855 passed every motion threshold
comfortably while carrying a coherent artifact affecting 46.4% of its white
matter orientations, against a cohort median of 19.1%. In a dataset whose
published QC covers only anatomical and functional data, an eigenvector-based
screen was the only barrier between that scan and the statistics. Given that
the finding here is a global FA shift of moderate size, a single badly
artifacted participant is not a negligible contribution.

**Pre-specification separates confirmation from search.** The corpus callosum
hypothesis was fixed in advance, so its failure is interpretable and the
analyses that followed are explicitly labelled exploratory. FDR correction
makes the nine regional results statistically defensible, but cannot make them
confirmatory: a finding located by searching 50 regions is a hypothesis
generated by the data, and the same data cannot also test it. Independent
replication is required before any of these regions should be treated as
established.

### 5.5 Future directions

**Replication.** The global finding needs testing in an independent cohort,
ideally acquired on a single scanner under a single protocol, which would
remove the confound that dominated cohort construction here and permit a larger
usable sample from the same source.

**Adequate power.** Detecting the corpus callosum effect observed here
(d = −0.29) with 80% power would require approximately 185 participants per
group. Future tests of focal hypotheses in this literature should be powered
for small effects rather than the moderate ones often assumed.

**Clinical variables.** Medication status, illness duration, age at onset and
symptom severity are available in ds000030 but were not modelled. Antipsychotic
exposure is a plausible contributor to white matter differences and cannot be
separated from illness effects in a cross-sectional design without it, although
ENIGMA detected no significant effect of medication dosage [3].

**Distinguishing detectability from pathology.** Because a global shift makes
the ranking of significant regions partly a map of measurement quality (§5.3),
a useful extension would be to model regional detectability explicitly — for
example by regressing each region's observed effect size against its skeleton
size, test–retest reliability and registration variability, and asking whether
any region departs from what detectability alone predicts. Regions exceeding
that baseline would be the genuine candidates for selective involvement.

**Acquisition permitting distortion correction.** A reverse-phase-encoded b = 0
acquisition costs under a minute of scanner time and would allow susceptibility
correction, making orbitofrontal, anterior temporal and brainstem regions
interpretable. The brainstem tracts ranking high in §4.3 are precisely those
this study cannot assess confidently.

**Models beyond the single tensor.** This study used one shell (64 directions
at b = 1000). Multi-shell acquisition, sampling several b-values, probes
different diffusion length scales and supports models the single tensor cannot.
Constrained spherical deconvolution [23] estimates a full orientation
distribution per voxel and so resolves crossing fibres, which otherwise lower
FA for purely geometric reasons. NODDI [24] separates intra-axonal,
extra-axonal and free-water compartments, distinguishing neurite density from
orientation dispersion — two microstructural properties that FA confounds.
Free-water elimination [25] models and removes a fast-diffusing extracellular
pool. This last is the most directly relevant: a diffuse FA reduction is
equally consistent with microstructural change and with a mild global increase
in extracellular water, and the present data cannot distinguish them.

**Constrained tensor estimation.** The degenerate fits in §3.4 arise from an
unconstrained least-squares solution. Positivity-constrained nonlinear fitting
would eliminate them at source rather than requiring voxel exclusion.

---

## 6. Limitations

### 6.1 Measurement

**Susceptibility distortion is uncorrected** (§2.5). Geometry is unreliable in
orbitofrontal, anterior temporal and brainstem regions. Because phase encoding
was uniform across participants this cannot create a spurious group difference,
but it means the brainstem tracts ranking high in §4.3 — pontine crossing
tract, cerebral peduncle — cannot be assessed confidently.

**FA is a proxy, not a myelin measure** (§1.3), and crossing fibres reduce it
for geometric reasons unrelated to tissue integrity. The same applies to AD and
RD, so the RD-dominant pattern in §4.4 is described rather than mechanistically
interpreted.

**Free water is unmodelled.** A diffuse FA reduction is equally consistent with
microstructural change and with a mild global increase in extracellular water.
The single-shell acquisition cannot separate them, which makes this the most
important ambiguity in the main finding.

**All diffusion QC is our own** (§2.9). The dataset's published QC covers
anatomical and functional data only.

### 6.2 Statistical power

With 36 participants per group, power to detect an effect of the magnitude
typically reported in this literature (d ≈ 0.4) is approximately 40%; 80% power
would require d ≈ 0.64. The global effect (d = −0.58) was therefore close to
adequately powered, while the corpus callosum test (d = −0.29) was not.

This is why a single primary hypothesis was pre-specified, why effect sizes and
confidence intervals accompany every p-value, and why the primary null is
presented as weak evidence of absence rather than as a demonstration that no
difference exists.

### 6.3 Sample

Participants excluded for protocol mismatch were disproportionately male
patients scanned on one scanner, so the final sample is not fully
representative of the original patient group. Clinical variables were available
but not modelled (§5.5).

**Residual ghosting imbalance.** The T1 ghosting flag remains imbalanced across
groups (29/7 control versus 21/15 schizophrenia, p = 0.072). Anatomical images
were excluded from the pipeline, so this cannot affect registration. However,
if the cause — possibly a headset worn during scanning — was present throughout
the session, the diffusion acquisition may also have been affected, and a
subtle global signal contamination imbalanced across groups could in principle
produce the global FA difference reported here.

A sensitivity analysis tested this directly. Global FA did not differ by
ghosting flag (0.5124 flagged versus 0.5111 unflagged), the flag did not
predict global FA when entered as an additional covariate (β = +0.0035,
p = 0.42), and the diagnosis effect was unchanged by its inclusion
(β = −0.0098, p = 0.013, versus β = −0.0090, p = 0.018 without it). The
ghosting imbalance therefore does not account for the global FA difference. The
pre-specified model was retained as primary and the flag was not added as a
routine covariate, since adding covariates on the basis that they strengthen a
result would undermine the pre-specification.

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
| `derivatives/roi_values.tsv` | Per-participant, per-region FA/MD/AD/RD |
| `derivatives/roi_stats.tsv` | Regional group statistics with FDR correction |
| `code/denoise_one.py` | MP-PCA denoising and Gibbs removal |
| `code/process_one_subject.sh` | Per-participant preprocessing |
| `code/run_batch.sh` | Parallel batch driver |
| `code/phase4f_motion_qc.py` | Motion QC and exclusion |
| `code/phase5_dtifit.sh` | Tensor fitting |
| `code/phase5_colorfa_check.py` | Orientation QC |
| `code/phase5b_crossvalidate.py` | Independent fit comparison |
| `validate_fit.py` | From-scratch tensor fit |
| `code/phase6a_register.sh` | Registration to FMRIB58_FA |
| `code/phase6b_roi_extract.py` | JHU atlas extraction and group models |
| `code/phase6c_interrogate.py` | Global/regional, quality and robustness checks |
| `code/phase8_tbss.sh` | TBSS pipeline, design matrix and permutation testing |
| `code/phase8b_tbss_report.py` | TBSS result summary and figures |
| `code/report_tables.py` | Report table generation |

All analysis subject lists are derived from `derivatives/motion_qc.tsv`
filtered on `exclude_final`, never from directory listings, since the
`raw/` and `derivatives/dti/` trees contain participants not in the analysis
cohort.

---

## References

1. Friston KJ, Frith CD. Schizophrenia: a disconnection syndrome? *Clinical
   Neuroscience* 1995;3(2):89–97. PMID: 7583624.
2. Poldrack RA, Congdon E, Triplett W, Gorgolewski KJ, Karlsgodt KH, Mumford
   JA, Sabb FW, Freimer NB, London ED, Cannon TD, Bilder RM. A phenome-wide
   examination of neural and cognitive function. *Scientific Data*
   2016;3:160110. doi:10.1038/sdata.2016.110.
3. Kelly S, Jahanshad N, Zalesky A, Kochunov P, et al. Widespread white matter
   microstructural differences in schizophrenia across 4322 individuals:
   results from the ENIGMA Schizophrenia DTI Working Group. *Molecular
   Psychiatry* 2018;23(5):1261–1269. doi:10.1038/mp.2017.170.
4. Friston K, Brown HR, Siemerkus J, Stephan KE. The dysconnection hypothesis
   (2016). *Schizophrenia Research* 2016;176(2–3):83–94.
   doi:10.1016/j.schres.2016.07.014.
5. Stejskal EO, Tanner JE. Spin diffusion measurements: spin echoes in the
   presence of a time-dependent field gradient. *Journal of Chemical Physics*
   1965;42(1):288–292. doi:10.1063/1.1695690.
6. Basser PJ, Mattiello J, LeBihan D. MR diffusion tensor spectroscopy and
   imaging. *Biophysical Journal* 1994;66(1):259–267.
   doi:10.1016/S0006-3495(94)80775-1.
7. Mori S, Wakana S, Nagae-Poetscher LM, van Zijl PCM. *MRI Atlas of Human
   White Matter*. Elsevier, Amsterdam, 2005.
8. Wakana S, Caprihan A, Panzenboeck MM, et al. Reproducibility of quantitative
   tractography methods applied to cerebral white matter. *NeuroImage*
   2007;36(3):630–644. doi:10.1016/j.neuroimage.2007.02.049.
9. Hua K, Zhang J, Wakana S, et al. Tract probability maps in stereotaxic
   spaces: analyses of white matter anatomy and tract-specific quantification.
   *NeuroImage* 2008;39(1):336–347. doi:10.1016/j.neuroimage.2007.07.053.
10. Smith SM, Jenkinson M, Johansen-Berg H, Rueckert D, Nichols TE, Mackay CE,
    Watkins KE, Ciccarelli O, Cader MZ, Matthews PM, Behrens TEJ. Tract-based
    spatial statistics: voxelwise analysis of multi-subject diffusion data.
    *NeuroImage* 2006;31(4):1487–1505. doi:10.1016/j.neuroimage.2006.02.024.
11. Smith SM, Jenkinson M, Woolrich MW, et al. Advances in functional and
    structural MR image analysis and implementation as FSL. *NeuroImage*
    2004;23(S1):S208–S219. doi:10.1016/j.neuroimage.2004.07.051.
12. Garyfallidis E, Brett M, Amirbekian B, et al. DIPY, a library for the
    analysis of diffusion MRI data. *Frontiers in Neuroinformatics* 2014;8:8.
    doi:10.3389/fninf.2014.00008.
13. Veraart J, Novikov DS, Christiaens D, Ades-aron B, Sijbers J, Fieremans E.
    Denoising of diffusion MRI using random matrix theory. *NeuroImage*
    2016;142:394–406. doi:10.1016/j.neuroimage.2016.08.016.
14. Kellner E, Dhital B, Kiselev VG, Reisert M. Gibbs-ringing artifact removal
    based on local subvoxel-shifts. *Magnetic Resonance in Medicine*
    2016;76(5):1574–1581. doi:10.1002/mrm.26054.
15. Smith SM. Fast robust automated brain extraction. *Human Brain Mapping*
    2002;17(3):143–155. doi:10.1002/hbm.10062.
16. Andersson JLR, Sotiropoulos SN. An integrated approach to correction for
    off-resonance effects and subject movement in diffusion MR imaging.
    *NeuroImage* 2016;125:1063–1078. doi:10.1016/j.neuroimage.2015.10.019.
17. Andersson JLR, Graham MS, Zsoldos E, Sotiropoulos SN. Incorporating outlier
    detection and replacement into a non-parametric framework for movement and
    distortion correction of diffusion MR images. *NeuroImage*
    2016;141:556–572. doi:10.1016/j.neuroimage.2016.06.058.
18. Jenkinson M, Bannister P, Brady M, Smith S. Improved optimization for the
    robust and accurate linear registration and motion correction of brain
    images. *NeuroImage* 2002;17(2):825–841. doi:10.1006/nimg.2002.1132.
19. Andersson JLR, Jenkinson M, Smith S. Non-linear registration, aka spatial
    normalisation. FMRIB Technical Report TR07JA2, University of Oxford, 2007.
20. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical
    and powerful approach to multiple testing. *Journal of the Royal
    Statistical Society, Series B* 1995;57(1):289–300.
21. Winkler AM, Ridgway GR, Webster MA, Smith SM, Nichols TE. Permutation
    inference for the general linear model. *NeuroImage* 2014;92:381–397.
    doi:10.1016/j.neuroimage.2014.01.060.
22. Smith SM, Nichols TE. Threshold-free cluster enhancement: addressing
    problems of smoothing, threshold dependence and localisation in cluster
    inference. *NeuroImage* 2009;44(1):83–98.
    doi:10.1016/j.neuroimage.2008.03.061.
23. Tournier J-D, Calamante F, Connelly A. Robust determination of the fibre
    orientation distribution in diffusion MRI: non-negativity constrained
    super-resolved spherical deconvolution. *NeuroImage* 2007;35(4):1459–1472.
    doi:10.1016/j.neuroimage.2007.02.016.
24. Zhang H, Schneider T, Wheeler-Kingshott CA, Alexander DC. NODDI: practical
    in vivo neurite orientation dispersion and density imaging of the human
    brain. *NeuroImage* 2012;61(4):1000–1016.
    doi:10.1016/j.neuroimage.2012.03.072.
25. Pasternak O, Sochen N, Gur Y, Intrator N, Assaf Y. Free water elimination
    and mapping from diffusion MRI. *Magnetic Resonance in Medicine*
    2009;62(3):717–730. doi:10.1002/mrm.22055.

---

## Appendix A — Decision log summary

| Decision | Evidence | Consequence |
|---|---|---|
| Restrict to SCZ vs control | Study scope | Bipolar and ADHD groups unused |
| Match on acquisition protocol | Scanner × diagnosis χ² p = 0.0001 | 9 patients excluded as unmatchable |
| Exclude sub-50029 | 32 volumes rather than 65 | — |
| Exclude T1 from pipeline | Ghosting 40% SCZ vs 15% control | FA registered directly to FMRIB58_FA |
| No susceptibility correction | No reverse-PE run, no fieldmap | Frontal, temporal and brainstem geometry unreliable |
| Motion thresholds fixed before FA inspection | — | 4 exclusions |
| Exclude sub-50005 | Sole participant on `tr7_sl50` after partner failed QC | Consistent with unmatchable rule |
| Exclude sub-10855 | L–R orientation 46.4% vs median 19.1%, threshold 29.2% | Motion was normal; only orientation QC caught it |
| Retain OLS over WLS | WLS increased degenerate fits (710/671/575 → 746/710/833) | FA > 1 exclusion rule retained |
| Exclude FA > 1 voxels from ROI means | 0.27% of voxels, $\sqrt{3/2}$ degenerate-fit signature | Applied to all regional averages |
| Unpaired analysis | Pair identity not retained | Age, sex, motion as covariates |
| Corpus callosum pre-specified as primary | ~40% power; one confirmatory test only | All else exploratory |
| Primary ROI not redefined after null | Genu alone d = −0.44 vs pooled d = −0.29 | Discussed in §5.2, not re-tested |
| FDR across 50 regions | Exploratory regional analysis | 9 regions at q < 0.05 |
| TBSS design order verified against input order | `randomise` matches rows by position | Check enforced before inference |
| Ghosting flag not added as routine covariate | Sensitivity analysis: no effect (p = 0.42) | Pre-specified model retained |