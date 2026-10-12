# Same Brain Age, Different Cognition

**Operationalising cognitive resilience as deviation from brain age**

Analysis code for a master's thesis at the Department of Psychology,
Carl von Ossietzky Universität Oldenburg.

Author: Inayet Dincer
Supervisors: Prof. Dr Andrea Hildebrandt, Dr Daniel Kristanto
Data: AABC / HCP-Aging Release 2 (*n* = 1,380, aged 36–90)

---

## The question

Two people can have brains that look equally old for their age and still differ
substantially in cognitive performance. This thesis treats that discrepancy as
**cognitive resilience**: cognitive performance conditional on brain status and
demographic background.

Brain status is measured as the **brain-age gap** — the difference between an age
predicted from neuroimaging and chronological age — estimated separately from
structural morphometry and from resting-state functional connectivity.

Three questions follow:

| | question | model |
|---|---|---|
| **H1** | Does the gap relate to cognition at all? | latent *g* ~ BAG + age + sex + education |
| **H2** | What is associated with the gap itself? | BAG ~ X + age + sex + education |
| **H3** | What is associated with resilience? | latent *g* ~ BAG + X + age + sex + education |

In the terminology of Bocancea et al. (2021), H2 measures **resistance** — the
mismatch between risk factors and brain integrity — and H3 measures
**resilience** — the mismatch between brain status and cognition.

---

## Repository contents

This repository contains **analysis code only**. No data and no results are
published here; see [Data availability](#data-availability).

```
code/
  brain_age/     brain-age estimation, Python
    01_build_connectivity.py     assemble 379x379 connectivity matrices
    02_build_structural.py       assemble 1,500 morphometric features
    03_build_phenotype.py        merges, recodes, exclusions (single point of truth)
    04_functional_bag.py         functional brain age, nested CV
    05_structural_bag.py         structural brain age, nested CV
    06_merge_bags.py             join both gaps to the phenotype
    07_sensitivity_bag.py        alternative correction, gap stability
    08_depth_diagnostic.py       grid-boundary diagnostic

  sem/           measurement and structural models, R / lavaan
    01_cfa_cognitive_bifactor.R  bifactor-(S-1) model of cognition
    02_cfa_psychosocial_factors.R  well-being and social connectedness
    03_impute.R                  multiple imputation, 50 datasets
    03b_convergence.R            imputation diagnostics
    04_h1_base.R                 H1
    05_h2_single.R               H2 primary, one exposure per model
    05b_h2_joint.R               H2 sensitivity, all exposures together
    06_h3_single.R               H3 primary
    06b_h3_joint.R               H3 sensitivity
    08_multiplicity.R            Bonferroni and Benjamini-Hochberg

  figures/       figures and tables for the thesis and presentation

docs/
  PIPELINE.md    what runs in what order, and why each step exists
```

Each script opens with a plain-English block stating where it sits in the run
order and which question it serves.

---

## Method summary

**Brain-age estimation.** Age was predicted with LightGBM under strictly nested,
family-grouped cross-validation: 10 repeated partitions × 10 outer folds × 5
inner folds. Every data-dependent preprocessing step — the tangent-space
reference, the scaler, the PCA basis, the head-size coefficients — was re-estimated
at each level of the hierarchy, on training participants only.

**Functional features.** Correlation matrices were converted back from Fisher-z,
checked for positive definiteness, and projected into the tangent space at the
Karcher mean of the training matrices. Correlation matrices are symmetric
positive-definite and so do not lie in a flat Euclidean space; the tangent-space
parameterisation respects that geometry and outperforms raw correlations in
benchmark comparisons (Dadi et al., 2019; Pervaiz et al., 2020). The 71,631
upper-triangle features were standardised and reduced to 100 principal components.

**Structural features.** 1,500 features: cortical thickness, surface area,
cortical volume and a T1w/T2w myelin proxy across 360 Glasser regions, plus 60
subcortical volumes. Intracranial volume was regressed out of the 780
size-dependent features using training-fold coefficients.

**Age-bias correction.** Brain-age predictions compress toward the sample mean,
so the raw gap correlates with chronological age. Calibration coefficients were
estimated by cross-fitting *inside each outer-training set* and applied unchanged
to the held-out fold, following More et al. (2023). No test participant
contributes to the correction applied to them.

**Measurement.** Latent general cognitive performance was estimated with a
bifactor-(S−1) model over seven indicators, with the general factor anchored to
two processing-speed tasks.

**Missing data.** Multiple imputation by chained equations (50 datasets, 50
iterations) for the observed predictors; full-information maximum likelihood for
the latent indicators; pooled by Rubin's rules.

**Multiplicity.** H1 is a family of two under Bonferroni (α = .025). H2 and H3 are
separate families of twelve under Benjamini-Hochberg (*q* = .05), because they ask
different questions.

---

## Reproducing the analysis

Access to AABC Release 2 is required. See [Data availability](#data-availability).

```bash
# 1. Python environment
pip install -r requirements.txt

# 2. Brain-age pipeline, in order
python code/brain_age/01_build_connectivity.py
python code/brain_age/02_build_structural.py
python code/brain_age/03_build_phenotype.py
python code/brain_age/04_functional_bag.py      # several hours
python code/brain_age/05_structural_bag.py      # several hours
python code/brain_age/06_merge_bags.py

# 3. Measurement models, imputation, hypothesis tests
Rscript code/sem/01_cfa_cognitive_bifactor.R
Rscript code/sem/02_cfa_psychosocial_factors.R
Rscript code/sem/03_impute.R                    # about 20 minutes
Rscript code/sem/03b_convergence.R
Rscript code/sem/04_h1_base.R
Rscript code/sem/05_h2_single.R
Rscript code/sem/06_h3_single.R
Rscript code/sem/08_multiplicity.R
```

### The run tag

**Every reported result comes from one run: `pca_it50`.** Use that tag and
nothing else to reproduce the thesis.

`FUNCTIONAL_TAG` is read by scripts in both languages, but it names two
different things depending on where you are in the chain, and this is the one
place the pipeline is easy to misread:

| where | what the variable names | values |
|---|---|---|
| `03_impute.R`, `03b_convergence.R` | the functional **pipeline** | `pca` (pre-specified) or `nopca` (sensitivity) |
| `04`–`08`, and the figure scripts | the **imputation run** | `pca_it50` |

`03_impute.R` combines the pipeline name with its iteration count and writes
its output under that combined name, so the file it produces for the
pre-specified pipeline is `mice_pca_it50_*.rds`. Everything downstream
addresses files by the **run** name. The R scripts from `04` onward therefore
default to `pca_it50`, not to `pca`.

```bash
# the reported analysis
Rscript code/sem/04_h1_base.R                      # defaults to pca_it50

# the unreduced-features sensitivity run
FUNCTIONAL_TAG=nopca python code/brain_age/06_merge_bags.py
```

**Why 50 iterations and not the 20 in the analysis plan.** At 20 iterations
perceived stress had a Gelman–Rubin PSRF of 1.07, above the 1.05 criterion; at
50 it was 1.02, with every other imputed variable below 1.01. Only the
iteration count changed — no method, no variable, no model. A 20-iteration run
exists in the project history and is superseded; it is not reproduced here and
no reported number comes from it.

See [`docs/PIPELINE.md`](docs/PIPELINE.md) for the full run order including
sensitivity analyses.

### Before you run anything

The R scripts set the project root as an absolute path near the top of each file:

```r
ROOT <- "C:/Users/inaye/Desktop/MASTERTHESIS"
```

Change this to your own path in each script under `code/sem/`, or set it once in
your `.Rprofile`. The Python scripts resolve their paths relative to the script
location and need no change.

The code is published exactly as it was run for the thesis, so this path is left
as it was rather than refactored after the fact.

---

## Requirements

**Python 3.11+** — numpy, pandas, scipy, scikit-learn, lightgbm, matplotlib.
Full versions in `requirements.txt`.

**R 4.3+** — lavaan, semTools, mice, dplyr.

---

## Data availability

AABC / HCP-Aging Release 2 is **registered-access** data, distributed by the
NIMH Data Archive under a data use agreement.

Neither the raw data nor any participant-level derived file may be redistributed.
Accordingly this repository contains:

- ✅ all analysis code
- ✅ documentation of every methodological decision
- ❌ no raw or processed data
- ❌ no participant-level results
- ❌ no imputation objects or model objects

The `.gitignore` enforces this. Researchers with their own approved access can
reproduce the analysis from the scripts here.

Access: <https://nda.nih.gov/>

---

## Key references

Bocancea, D. I., van Loenhoud, A. C., Groot, C., Barkhof, F., van der Flier, W. M.,
& Ossenkoppele, R. (2021). Measuring resilience and resistance in aging and
Alzheimer disease using residual methods. *Neurology, 97*(10), 474–488.

Dadi, K., Rahim, M., Abraham, A., Chyzhyk, D., Milham, M., Thirion, B., &
Varoquaux, G. (2019). Benchmarking functional connectome-based predictive models
for resting-state fMRI. *NeuroImage, 192*, 115–134.

More, S., Antonopoulos, G., Hoffstaedter, F., Caspers, J., Eickhoff, S. B., &
Patil, K. R. (2023). Brain-age prediction: A systematic comparison of machine
learning workflows. *NeuroImage, 270*, 119947.

Pervaiz, U., Vidaurre, D., Woolrich, M. W., & Smith, S. M. (2020). Optimising
network modelling methods for fMRI. *NeuroImage, 211*, 116604.

Stern, Y., Gazes, Y., Razlighi, Q., Steffener, J., & Habeck, C. (2017). Cognitive
reserve and brain maintenance: Orthogonal concepts in theory and practice.
*Cerebral Cortex, 27*(8), 3962–3969.

---

## Licence

Code released under the MIT Licence. See [`LICENSE`](LICENSE).

The licence covers this code only. It confers no rights to AABC / HCP-Aging data.
