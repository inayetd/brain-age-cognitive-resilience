"""
STEP 6, Merge the two brain age gaps with the phenotype into one analysis table.

Writes data/clean/analysis_table_<TAG>.csv, the single file the R / lavaan models read.
Also prints the H4a descriptive quantities: the correlation between the two gaps, their
agreement, and the size of the intersection sample.

INPUTS   bag_functional_pca.csv (or _nopca), bag_structural.csv
         ../../data/clean/analysis_phenotype.csv
OUTPUT   ../../data/clean/analysis_table_<TAG>.csv

FUNCTIONAL_TAG selects which functional run is merged. 'pca' is the primary specification;
'nopca' is the pre-specified sensitivity run on the unreduced tangent features. Set the
environment variable FUNCTIONAL_TAG=nopca to build the sensitivity table instead; the
tag is carried into the output filename, so the two runs cannot overwrite each other.
"""
# ============================================================================
# Step 6 of 8. Joins the two gaps to the phenotype.
# Produces the one table every R script reads. Also reports how far the two gaps
# agree with each other, which is the descriptive part of H4.
# Run before: any script in code/sem/.
# ============================================================================
import os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr

N_BOOT = 2000   # family-clustered bootstrap replicates for the H4a interval

WORK = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(WORK, '..', '..'))
OUTD = os.path.join(ROOT, 'data', 'clean')

ph = pd.read_csv(os.path.join(OUTD, 'analysis_phenotype.csv'))
FUNCTIONAL_TAG = os.environ.get('FUNCTIONAL_TAG', 'pca')
fn = os.path.join(WORK, f'bag_functional_{FUNCTIONAL_TAG}.csv')
print(f'functional run merged: {FUNCTIONAL_TAG}')
sn = os.path.join(WORK, 'bag_structural.csv')

out = ph.copy()
for path, tag in ((fn, 'functional'), (sn, 'structural')):
    if not os.path.exists(path):
        print(f'  {tag} BAG not found ({os.path.basename(path)}), skipping')
        continue
    cols = ['participant_id', f'BAG_{tag}', 'BAG_sd_across_seeds']
    src = pd.read_csv(path)
    if 'BAG_n_valid_repetitions' in src.columns:
        cols.append('BAG_n_valid_repetitions')
    b = src[cols].rename(columns={'BAG_sd_across_seeds': f'BAG_{tag}_sd',
                                  'BAG_n_valid_repetitions': f'BAG_{tag}_n_valid'})
    out = out.merge(b, on='participant_id', how='left')

hasF = 'BAG_functional' in out.columns
hasS = 'BAG_structural' in out.columns
if hasF and hasS:
    both = out[['BAG_functional', 'BAG_structural']].notna().all(axis=1)
    out['both_modalities'] = both
    a = out.loc[both, 'BAG_functional'].to_numpy()
    b = out.loc[both, 'BAG_structural'].to_numpy()
    # Lin's concordance correlation coefficient (CCC), not an ICC, though the two
    # coincide under some assumptions.
    ccc = 2 * np.cov(a, b)[0, 1] / (a.var() + b.var() + (a.mean() - b.mean()) ** 2)
    r = np.corrcoef(a, b)[0, 1]
    n = len(a)

    # Family-clustered bootstrap interval rather than the Fisher z formula, which
    # assumes independent observations. Resampling families keeps relatives together,
    # matching the fold structure.
    fam = out.loc[both, 'pedigree'].astype(str).to_numpy() if 'pedigree' in out.columns else \
          np.array([str(i) for i in range(len(a))])
    if 'pedigree' not in out.columns:
        sys.stderr.write('WARNING: no pedigree column; the H4a bootstrap treats every '
                         'participant as unrelated, which will understate the interval.\n')
    fams = np.unique(fam)
    idx_by_fam = {f: np.where(fam == f)[0] for f in fams}
    rng = np.random.default_rng(0)
    boot = []
    for _ in range(N_BOOT):
        pick = rng.choice(fams, size=len(fams), replace=True)
        ix = np.concatenate([idx_by_fam[f] for f in pick])
        if len(ix) > 2:
            rb = np.corrcoef(a[ix], b[ix])[0, 1]
            if np.isfinite(rb):
                boot.append(rb)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    print(f'\nH4a descriptives (intersection sample n={int(both.sum())}, '
          f'{len(fams)} families):')
    print(f'  Pearson r                 {r:.3f}   95% CI [{lo:.3f}, {hi:.3f}]  '
          f'(family-clustered bootstrap, {len(boot)} of {N_BOOT} replicates)')
    # H4 is exploratory, so no support/rejection verdict is printed.
    print(f'  Spearman rho              {spearmanr(a, b).statistic:.3f}')
    print(f"  Lin's CCC                 {ccc:.3f}")
    print(f'  shared variance           {r ** 2:.3f}')
    print(f'  sign agreement            {100 * np.mean(np.sign(a) == np.sign(b)):.1f}%  '
          f'(descriptive only, sign is unstable near zero)')
    # Deciles of each gap, to report tail overlap
    for q, name in ((0.9, 'top'), (0.1, 'bottom')):
        fa = a >= np.quantile(a, q) if q > 0.5 else a <= np.quantile(a, q)
        fb = b >= np.quantile(b, q) if q > 0.5 else b <= np.quantile(b, q)
        print(f'  {name} decile overlap       {100 * (fa & fb).sum() / fa.sum():.1f}%')

# Separate tables per run: a single filename would let whichever ran last become the
# analysis table.
p = os.path.join(OUTD, f'analysis_table_{FUNCTIONAL_TAG}.csv')
out.to_csv(p, index=False)
print(f'\nwrote {p}   rows {len(out)}   columns {out.shape[1]}')
print(f'  SEM-eligible {int(out.sem_eligible.sum())}'
      + (f'   with both modalities {int(out.both_modalities.sum())}' if hasF and hasS else ''))
