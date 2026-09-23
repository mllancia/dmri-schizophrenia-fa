
---

## Phase 4f — Motion QC (2026-09-14)

Thresholds fixed **before any FA value was inspected**:

- mean volume-to-volume RMS displacement > 2.0 mm
- maximum volume-to-volume RMS displacement > 3.0 mm
- eddy-flagged outlier slices > 10%

Metrics extracted from `eddy_corrected.eddy_movement_rms` (column 2, relative
displacement) and `eddy_corrected.eddy_outlier_map`. Written to
`derivatives/motion_qc.tsv`.

**Group comparison at n = 78:** CONTROL 0.243 ± 0.142 mm, SCZ 0.267 ± 0.165 mm;
Welch t(75) = −0.70, p = 0.485; Mann–Whitney p = 0.589; Cohen's d = −0.16,
95% CI [−0.60, 0.29]. No significant difference. Motion retained as a GLM
covariate regardless: at this sample size the CI is consistent with a moderate
effect, and reviewers will expect it.

**Note on `max_rel_rms`.** Group means differ substantially (0.793 vs 1.152) but
medians barely do (0.560 vs 0.605), and the SCZ SD is more than double that of
controls. The difference is driven by individual extreme values, chiefly
sub-50058 at 9.04 mm, not by a group tendency. Report medians and SDs alongside
means.

**Excluded (4), all on the maximum-displacement criterion:**

| Subject | Group | mean_rel | max_rel | %outlier |
|---|---|---|---|---|
| sub-10193 | CONTROL | 0.69 | 3.55 | 1.3 |
| sub-11062 | CONTROL | 0.82 | 3.02 | 2.3 |
| sub-50022 | SCZ | 0.94 | 5.79 | 3.8 |
| sub-50058 | SCZ | 0.65 | 9.04 | 0.8 |

Cohort 78 → 74.

---

## Phase 4f — Protocol exclusion (2026-09-14)

**sub-50005 (SCZ) excluded.** Its matched control sub-10193 failed motion QC,
leaving sub-50005 the sole subject on acquisition protocol
`35343_tr7_sl50` with no control counterpart anywhere in the cohort.

Consistent with the Phase 2 rule that removed nine male SCZ subjects on scanner
35426: subjects with no protocol-matched counterpart in the opposite group are
excluded. Applying the rule in Phase 2 and abandoning it here would be
indefensible.

Cohort 74 → 73.

---

## Phase 5 — Orientation QC (2026-09-14)

Colour-FA maps (|V1| × FA) generated for all 74 fitted subjects and inspected
as a montage, then individually for any subject appearing atypical. Criterion:
corpus callosum red (L–R), corticospinal tract blue (S–I), SLF green (A–P).

**Cohort-wide result: orientation check passed.** No subject showed scrambled
or swapped colour assignments, confirming the rotated bvecs were applied
correctly throughout.

Five subjects were inspected individually at three slice levels. Four
(sub-10460, sub-10530, sub-50007, sub-50032) proved normal; their montage
appearance was explained by slice-level variation, since the montage samples
50% of the z dimension and subjects differ in slice count and head position.
**A montage panel is adequate for detecting gross orientation failure and
inadequate for judging anything subtler.**

**sub-10855 (CONTROL) excluded — coherent orientation artifact.**

Visual: broad red regions through frontal and occipital white matter,
indicating an extensive sheet of apparently coherent left–right diffusion with
no anatomical counterpart.

Quantified as the proportion of WM voxels (FA > 0.25) with |V1ₓ| > 0.8, across
all 74 subjects; outlier rule median + 3 × IQR fixed before inspecting the
ranking.

- Cohort median 19.1%, IQR 3.4% → threshold 29.2%
- sub-10855: **46.4%**, the only subject exceeding threshold
- Next highest: sub-10487 at 27.1%
- Mean diffusivity also elevated: 0.00111 vs cohort typical ≈ 0.00093 mm²/s
- Eddy outlier report: slice 1 flagged repeatedly across volumes at ≈ −4 SD,
  a systematic per-slice pattern rather than random motion

**Motion was entirely normal**: mean relative 0.304 mm, max 0.966 mm, 2.4%
outlier slices — comfortably within every threshold.

Metric showed no association with diagnosis (CONTROL 19.62 ± 5.28, SCZ
19.23 ± 2.64; the elevated CONTROL SD is attributable entirely to this
subject). Written to `derivatives/orientation_qc.tsv`.

**Methodological point for the report:** motion QC alone was insufficient. In a
dataset whose published QC (MRIQC) covers only anatomical and functional data,
an orientation-based screen was the only thing preventing a badly artifacted
scan from entering the statistics.

Cohort 73 → **72 (36 CONTROL, 36 SCZ)**.

---

## Phase 5 — Balance after exclusions (2026-09-14)

Re-verified at n = 72:

| Variable | CONTROL | SCZ | p |
|---|---|---|---|
| Scanner 35343 / 35426 | 22 / 14 | 22 / 14 | 1.000 |
| Protocol tr8.4_sl60 / tr9_sl60 | 22 / 14 | 22 / 14 | 1.000 |
| Sex F / M | 10 / 26 | 10 / 26 | 1.000 |
| Age, mean years | 34.2 | 35.8 | 0.440 |
| T1 ghost absent / present | 29 / 7 | 21 / 15 | 0.072 |

Scanner, protocol and sex are now **exactly** balanced. The protocol cell
`35343_tr7_sl50` no longer exists in the cohort.

T1 ghosting remains imbalanced but is moot: anatomical images were never
downloaded and registration goes FA → FMRIB58_FA directly (Route B). Reported
in the limitations as a flagged-and-avoided confound.

---

## Phase 5 — Degenerate tensor fits (2026-09-14)

**Rule: voxels with FA > 1 are excluded from all ROI averages.** Fixed before
any regional FA value was computed.

Maximum FA in every subject is exactly 1.224744 = √(3/2). This arises when
fitted eigenvalues sum to zero, which happens when unconstrained OLS returns
negative eigenvalues — mathematically permitted, physically impossible. FSL's
FA formula reduces to exactly √(3/2) in that case.

- Affected voxels: 0.27% of in-mask voxels on average (range 0.199–0.364%)
- ~2% of in-mask voxels have a negative smallest eigenvalue — low and normal
  for unconstrained linear fitting at b = 1000
- **Zero** fall in CSF-range MD; ~36% fall in WM-range MD (< 0.0012 mm²/s)
- ~40% survive eroding the brain mask by 3 voxels, so they are not purely a
  rim artifact and can fall inside atlas regions — which is what makes the
  exclusion rule load-bearing rather than cosmetic

**WLS sensitivity check.** `dtifit --wls` on three subjects left the maximum FA
unchanged to six decimal places (1.224744 / 1.224743 vs 1.224745). OLS
retained. `[VERIFY]` see the voxel-count comparison in
`derivatives/wls_comparison.txt` — the maximum is insensitive to how many
voxels were fixed, so the report should cite the count, not the max.

---

## Analysis design notes

- `cohort_final.tsv` has no `pair_id` column; pair identity was never retained.
  Matching served to balance protocol, sex and age at the **group** level.
  All analyses are unpaired by design, with age, sex and mean relative
  displacement as covariates.
- **All subject lists must be built from `derivatives/motion_qc.tsv` filtered
  on `exclude_final`.** Never from `ls derivatives/dti/sub-*` or `raw/sub-*`:
  `raw/` holds ~110 subjects, `derivatives/dti/` holds 74, the analysis cohort
  is 72.
- Corpus callosum pre-specified as the single primary ROI hypothesis.
  Everything else is exploratory. Power ≈ 40% at d = 0.4; d ≈ 0.64 needed for
  80%.

**WLS sensitivity result (2026-09-14).** `dtifit --wls` did NOT reduce
degenerate fits; voxel counts with FA > 1 rose in all three test subjects
(OLS 710 / 671 / 575 -> WLS 746 / 710 / 833). Likely cause: FSL derives WLS
weights from an initial OLS fit, so in noise-dominated voxels the weights are
themselves corrupted and amplify the bad measurements. OLS retained.
Counts in `derivatives/wls_comparison.txt`. Note the earlier max-FA comparison
was uninformative - a single degenerate voxel pins the max at sqrt(3/2)
regardless of how many others changed.
