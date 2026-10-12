# The pipeline, step by step

What runs, in what order, why it exists, and what it produces. Read this before
reading any script.

The design principle throughout: **every quantity estimated from data is
estimated on training participants only, at every level of the hierarchy.** Most
of the complexity in these scripts exists to enforce that.

---

## Stage 1 — Build the inputs

### `01_build_connectivity.py`

The release ships resting-state connectivity one seed region at a time, as
Fisher-z values with the diagonal set to zero. This script reassembles them into
one 379 × 379 matrix per participant and records the region order and session
identifiers.

Nothing is modelled, transformed, or excluded here. Undoing the Fisher-z
transform happens later, in step 4, because that is where the matrix geometry
starts to matter.

*Produces:* `corr_V1.npy`, `meta.json`

### `02_build_structural.py`

Assembles 1,500 structural features from the released morphometry tables:
cortical thickness, surface area, cortical volume and a T1w/T2w myelin proxy
across 360 Glasser regions, plus 60 subcortical volumes. Also extracts
intracranial volume, needed in step 5.

The script records which 780 features are size-dependent (`volmask.npy`). Surface
area and volumes scale with head size; thickness and the myelin ratio do not.

*Produces:* `Xstruct.npy`, `icv.npy`, `volmask.npy`, feature names

### `03_build_phenotype.py`

**The single point where participants are selected.** Every merge, recode and
exclusion happens here: the age floor, the neurological exclusions, the
reverse-scoring of loneliness, the sign of Trail Making B.

This matters more than it looks. If the measurement model were fitted on one
sample and the hypothesis tests on another, the factor scores would not describe
the people they are used on. Concentrating every exclusion in one script makes
that impossible by construction, and downstream scripts never reopen a raw
release file.

*Produces:* `data/clean/analysis_phenotype.csv`

---

## Stage 2 — Estimate brain age

Steps 4 and 5 share an identical cross-validation structure, deliberately. The
two modalities must differ only in their features, or any comparison between them
would be confounded by the pipeline.

### The cross-validation hierarchy

```
10 repeated partitions
└── 10 outer folds              held out; never seen during fitting
    ├── 5 inner folds           hyperparameter selection only
    └── 5 calibration folds     age-bias coefficients only
```

**Outer folds** are cut over *families*, not participants, and stratified on
family median age and modal sex. Only 16 of 1,368 families have more than one
member, so this changes little — but relatives have correlated brain structure,
and a sibling pair straddling a fold boundary means the model is partly
memorising.

**Inner folds** are the only place hyperparameters are chosen: learning rate
(.03, .05, .10) crossed with maximum tree depth (3, 5, 7). The selected
configuration is then refit on the full outer-training set.

**Preprocessing is refitted at every level.** In both scripts this is enforced by
a single function, `fit_preprocessing()`, which is called three times per outer
fold with three different fit sets. If you are checking this code for leakage,
that function and its call sites are where to look.

### `04_functional_bag.py` — the functional gap

```
Fisher-z values
  → tanh, restore unit diagonal         recover the correlation matrix
  → positive-definiteness check         shrinkage toward identity if it fails
  → tangent space at the Karcher mean   estimated from training participants
  → upper triangle, 71,631 features
  → standardise                         training mean and SD
  → PCA, 100 components                 training basis
  → LightGBM
  → age-bias correction
  → BAG
```

**Why tangent space.** A correlation matrix is symmetric positive-definite, and
those do not form a flat vector space. Averaging two of them the ordinary way can
produce a misleading answer: a matrix with *r* = +.9 and one with *r* = −.9
average to the identity — zero correlation, more apparent independence than
either input had. The tangent-space projection expresses each participant
relative to a geometry-aware centre instead. Benchmark comparisons across many
cohorts and pipelines favour it over raw or partial correlations (Dadi et al.,
2019; Pervaiz et al., 2020).

**The Karcher mean** is that centre. It is found by iteration: start from the
ordinary average, express every participant as a deviation from that guess using
the matrix logarithm, average the deviations, update the guess with the matrix
exponential, repeat. At the true centre the average deviation is zero, which is
also the convergence criterion the script logs. For single numbers this procedure
reduces exactly to the geometric mean.

**The reference is re-estimated per fold**, like everything else.

*Produces:* `bag_functional_<TAG>.csv`, `bag_metrics_<TAG>.json`

### `05_structural_bag.py` — the structural gap

Same machinery, three differences:

- no tangent projection — the features are already a plain matrix
- no PCA — 1,500 features is manageable for a tree ensemble
- **head-size adjustment** — intracranial volume regressed out of the 780
  size-dependent features, with coefficients from training participants only.
  Without it the model can learn head size, which has nothing to do with ageing.

*Produces:* `bag_structural.csv`, `bag_structural_metrics.json`

### The age-bias correction, in both scripts

Brain-age predictions regress toward the sample mean: young participants are
predicted too old, old participants too young. The raw gap therefore correlates
negatively with chronological age, and without correction every downstream
result would be partly an age effect in disguise.

The naive fix — fitting the correction on the predictions being corrected — is
circular. Instead:

```python
# inside the outer-training set only
for calibration_fold in GroupKFold(5).split(train):
    fit on 4/5, predict the held-out 1/5      # out-of-sample training predictions
beta, alpha = polyfit(age_train, those_predictions)

corrected = (prediction - alpha) / beta       # applied to the test fold
BAG       = corrected - chronological_age
```

α and β are estimated without any test participant, then applied unchanged. This
follows More et al. (2023).

One consequence worth knowing: correction **reduces** apparent accuracy, because
dividing by β < 1 rescales the predictions outward and increases their variance.
Accuracy is therefore reported from the *uncorrected* predictions, which is the
honest number — correcting with a transformation derived from age and then
reporting accuracy against age would be partly circular.

A guard refuses to correct any fold whose β falls to 0.2 or below, since dividing
by a near-zero slope would turn small prediction differences into enormous gaps.
In the completed runs this guard never fired.

### `06_merge_bags.py`

Joins both gaps to the phenotype and produces the one table every R script reads.
Also reports how far the two gaps agree — the descriptive component of H4.

*Produces:* `data/clean/analysis_table_<TAG>.csv`

---

## Stage 3 — Measurement models

### `01_cfa_cognitive_bifactor.R`

A bifactor-(S−1) model over seven indicators:

```
g    =~ flanker + card_sort + list_sorting + pattern_comp + picture_seq + ravlt + trails_b_rev
Exec =~ a*flanker + a*card_sort
Mem  =~ picture_seq + ravlt + list_sorting
g ~~ 0*Exec    g ~~ 0*Mem    Exec ~~ 0*Mem
```

**Why S−1 rather than a symmetric bifactor.** Symmetric bifactor models are prone
to anomalous results — vanishing or negative specific-factor variances, and a
general factor whose meaning shifts with which indicators happen to be included.
The S−1 form anchors the general factor to a reference domain. Here pattern
comparison and Trails-B load on *g* only, making it processing-speed-anchored, so
the two specific factors describe what executive and memory tasks share *beyond*
speed — an interpretable quantity rather than a floating residual.

The two executive loadings are constrained equal because a two-indicator factor
is not otherwise identified.

### `02_cfa_psychosocial_factors.R`

Two correlated factors: psychological well-being (three indicators) and social
connectedness (four, with loneliness reverse-scored). Perceived stress has a
single indicator and so enters the structural models as an observed variable.

---

## Stage 4 — Missing data

### `03_impute.R`

Household income is 25% missing and carries the main finding, so this step
matters more than usual.

**Two mechanisms, deliberately:**

- **Chained equations** for the observed predictors — income, sleep quality,
  perceived stress, education, physical activity. 50 datasets, 50 iterations.
  Predictive mean matching for continuous variables, ordinal models for ordered
  ones.
- **Full-information maximum likelihood**, inside lavaan, for the latent
  indicators.

The cognitive and psychosocial indicators, and the brain-age gap, are included in
the imputation model as *auxiliary* variables and then **restored to missing**
afterwards. This satisfies congeniality — a variable in the analysis model must be
in the imputation model — without any imputed outcome ever entering an analysis.

The 50 iterations (rather than the default 5) followed a convergence warning on
an earlier run. `03b_convergence.R` reports the diagnostics.

Everything downstream pools across all 50 datasets by Rubin's rules:
*T* = *W̄* + (1 + 1/*m*)*B*.

---

## Stage 5 — Hypothesis tests

All models: `estimator = "MLR"`, `missing = "fiml"`, `std.lv = TRUE`.
Standardisation is applied *within* each imputed dataset, not on pooled data.

### `04_h1_base.R` — does the gap relate to cognition?

```
g    ~ BAG + age + sex + education
Exec ~ BAG + age + sex + education
Mem  ~ BAG + age + sex + education
```

The specific factors are regressed on the same predictors so that the
*specificity* of the association is observed rather than assumed. If the gap had
related to memory specifically rather than to general ability, that would be a
substantively different finding — and a model containing only *g* could not have
detected it.

The gap is standardised here.

### `05_h2_single.R` — what relates to the gap itself?

One model per exposure. Six exposures × two modalities = twelve models.

```
BAG ~ X + age + sex + education
```

**Why one at a time.** H2 asks whether each exposure is associated with the gap
and how large that association is. Entering all six simultaneously answers a
different question — what is unique to each once the others are partialled out —
and with correlated exposures the joint coefficients can be smaller than any
bivariate association without any of them being null.

The gap stays **in years** here, so coefficients read directly as years of brain
age per standard deviation of the exposure.

**Which specification is confirmatory.** `05b_h2_joint.R` enters all six
exposures simultaneously, and **that is the confirmatory test**: the
preregistration (Dincer, Kristanto & Hildebrandt, September 2025) states that
"the predictor coefficients from the simultaneous H2 and H3 models will
constitute the confirmatory tests" and that single-predictor models "will be
reported only as supplementary descriptive analyses". `05_h2_single.R` is
therefore the supplementary one, despite what its own header says — see the
known deviation noted under Stage 6.

### `06_h3_single.R` — what relates to resilience?

```
g ~ BAG + X + age + sex + education
```

By the **Frisch–Waugh–Lovell theorem**, the coefficient on X in this model is
exactly its association with the part of cognition that the gap and the
demographics do not explain — which is what this thesis means by resilience.

**No separate residual score is computed or used.** Computing a residual and then
regressing it on predictors is a two-step procedure that ignores first-step
uncertainty and produces standard errors that are too small. The one-step
conditional model gives the same point estimate with correct inference, and lets
the measurement model carry the reliability of the outcome.

The gap is standardised here, unlike in H2.

`06b_h3_joint.R` enters all six exposures simultaneously and, as for H2, **is
the confirmatory specification**. `06_h3_single.R` is supplementary and
descriptive.

---

## Stage 6 — Multiplicity

### `08_multiplicity.R`

Reads the three estimate files, logs the MD5 checksum of each, and applies the
pre-specified corrections.

| family | tests | correction |
|---|---|---|
| H1 | 2 | Bonferroni, α = .025 |
| H2 | 12 | Benjamini-Hochberg, *q* = .05 |
| H3 | 12 | Benjamini-Hochberg, *q* = .05 |

**Why the corrections differ.** H1 is a confirmatory test of one hypothesis in two
modalities; a false positive in either would be claimed as support, so the family
is small and family-wise control is appropriate — and over two tests it costs
almost nothing. H2 and H3 are screening questions over six exposures with no
strong prior about which. Bonferroni over twelve would control any false positive
at the cost of missing most true ones, which is the wrong trade for a shortlist;
Benjamini-Hochberg controls the expected *proportion* of false discoveries among
the rejections.

**Why H2 and H3 are separate families.** They ask different questions. Pooling
them would mean a finding about cognition changes its status according to how
many lifestyle variables happened to be measured. Family definitions should
follow the inferential structure, not the spreadsheet.

### Known deviation — which rows `08_multiplicity.R` corrects

`08_multiplicity.R` as committed reads `h2_single_estimates_*.csv` and
`h3_estimates_*.csv`, that is the **single-predictor** models, and the headers of
`05b_h2_joint.R` and `06b_h3_joint.R` describe themselves as sensitivity
analyses. Both reflect an earlier draft of the analysis plan, in which the
single-predictor models were confirmatory.

The registered plan that was followed reverses this: the twelve tests in each
Benjamini-Hochberg family are the **simultaneous** coefficients, and the
single-predictor models are supplementary. The confirmatory *q* values reported
in the thesis are therefore computed over the twelve joint coefficients, not
over the rows in `multiplicity_*.csv`.

This changes no conclusion. Household income is the only exposure to survive
correction, in both modalities and in both outcomes, under either family
definition:

| | single | joint |
|---|---|---|
| H2 structural | −0.98 (*q* < .001) | −0.93 (*q* < .001) |
| H2 functional | −0.82 (*q* = .038) | −0.85 (*q* = .028) |
| H3 structural | +0.104 (*q* = .004) | +0.099 (*q* = .008) |
| H3 functional | +0.119 (*q* = .001) | +0.114 (*q* = .003) |

The scripts are left byte-identical to the versions that produced the results,
because each one records its own MD5 in its provenance file and editing a
comment would break that chain. The deviation is recorded here instead.

---

## Sensitivity analyses

| script | what it varies |
|---|---|
| `07_sensitivity_bag.py` | residualised instead of rescaled bias correction; per-participant gap stability |
| `08_depth_diagnostic.py` | whether the tree-depth grid was wide enough |
| `05b_h2_joint.R` | all exposures entered simultaneously |
| `06b_h3_joint.R` | same, for resilience |
| `FUNCTIONAL_TAG=nopca` | the full 71,631 tangent features instead of 100 components |

---

## Conventions

**The run tag.** Every reported result comes from one run, `pca_it50`.

`FUNCTIONAL_TAG` is read from the environment in both languages, but it names
two different things depending on where you are in the chain:

| where | what it names | values |
|---|---|---|
| `03_impute.R`, `03b_convergence.R` | the functional **pipeline** | `pca` (pre-specified) or `nopca` (sensitivity) |
| `04`–`08`, figure scripts | the **imputation run** | `pca_it50` |

`03_impute.R` appends the iteration count to the pipeline name and writes
`mice_pca_it50_*.rds`. Everything downstream addresses files by that combined
run name, which is why the R scripts from `04` onward default to `pca_it50`
rather than `pca`. Passing `FUNCTIONAL_TAG=pca` to `04` would send it looking
for files that do not exist.

The analysis plan specified 20 iterations. Perceived stress had not converged
there (PSRF 1.07 against a 1.05 criterion) and had at 50 (1.02), so the
reported run uses 50. Only the iteration count differs. The 20-iteration run is
superseded and no reported number comes from it.

**Provenance files.** Every R script writes a `*_provenance_*.txt` recording the
seed, package versions, and the checksum of the script itself. The multiplicity
script additionally records the checksum of each input, so any corrected value
can be traced to the exact result file it came from.

**Guards.** The R scripts contain `stop()` calls that halt the run if an
assumption fails — a gap with no missing values where missingness is expected,
a measurement model fitted on the wrong sample, an inadmissible solution. These
are deliberate. Do not submit these scripts to the console line by line, because
that lets execution continue past a guard that was meant to stop it.
