"""
STEP 8, Diagnostic only: was the tree-depth grid wide enough?

WHY THIS EXISTS
  The analysis plan fixes the inner-loop search to learning rate (.03, .05, .10)
  crossed with maximum depth (3, 5, 7). Across the three completed runs the inner
  loop chose depth 3 in every single fold:

      functional PCA     depth 3 in 100 of 100 folds
      functional no-PCA  depth 3 in 100 of 100 folds
      structural         depth 3 in 100 of 100 folds

  The learning rate did vary across all three of its values, so the selection is
  discriminating between settings rather than failing. It is depth specifically
  that is pinned, and it is pinned at the LOWEST value in the grid. When an
  optimum sits on the boundary of a search, the honest position is that you do
  not know whether something outside the boundary would have been better.

  This script asks that question by extending the grid downwards to depths 1 and
  2 and comparing inner-fold performance.

WHAT THIS IS NOT
  It does not change the primary analysis and must not be used to. The grid is
  pre-specified; widening it in the models that produce the reported brain age
  gaps would be a departure from the analysis plan. The result is reported as a
  post-hoc check:

      if depth 3 still wins, the boundary is real and that can be stated;
      if 1 or 2 wins, the models are slightly over-complex, and that is a
      limitation to report rather than a reason to refit.

  Nothing here is written to any file the primary pipeline reads.

HOW IT WORKS
  The split, the preprocessing and the scoring are imported from
  05_structural_bag.py rather than copied, so this diagnostic cannot drift away
  from the procedure it is diagnosing. Only the depth grid differs.

  The structural arm is used because it is the cheap one: no tangent-space
  projection, no Karcher mean, no PCA, so a fold costs seconds rather than
  minutes. Depth was pinned identically in both modalities, so the question is
  the same one.

RUNTIME
  Roughly 15 LightGBM fits per inner fold, 5 inner folds per outer fold. On this
  feature set that is a couple of minutes per outer fold. N_FOLDS below controls
  how many outer folds of seed 0 are examined; three is enough to see whether
  the ordering is stable, ten removes any doubt.

OUTPUT  ../../results/depth_diagnostic.csv
        and a table printed to the console
"""
# ============================================================================
# Step 8 of 8. Diagnostic only.
# The inner loop selected the smallest tree depth on the grid in every one of the
# 300 folds. This script asks whether an even smaller depth would have been chosen
# if it had been offered, which is the honest way to report a grid boundary.
# ============================================================================
import importlib.util
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
RESD = os.path.join(ROOT, 'results')

# Depths to compare. The pre-specified grid is (3, 5, 7); 1 and 2 are the values
# the search never reached and are the reason this script exists.
DEPTH_GRID = (1, 2, 3, 5, 7)
N_FOLDS    = 3        # outer folds of seed 0 to examine; raise to 10 for all
SEED       = 0


def load_primary():
    """Import 05_structural_bag.py without running it.

    The module name starts with a digit, so the normal import statement cannot
    reach it. Its main() is behind the usual __main__ guard, so importing is
    safe: the pipeline does not start.
    """
    path = os.path.join(HERE, '05_structural_bag.py')
    if not os.path.exists(path):
        raise SystemExit(f'not found: {path}')
    spec = importlib.util.spec_from_file_location('bag05', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['bag05'] = mod
    spec.loader.exec_module(mod)
    return mod


def load_data(m):
    """The same participants, features and exclusions as the primary run."""
    ids = [re.sub(r'_V\d+.*$', '', str(s))
           for s in json.load(open(os.path.join(HERE, 'meta.json')))['ids']]
    X   = np.load(os.path.join(HERE, 'Xstruct.npy'))
    icv = np.load(os.path.join(HERE, 'icv.npy'))
    vm  = np.load(os.path.join(HERE, 'volmask.npy'))

    ph = pd.read_csv(m.PHENO).set_index('participant_id')
    keep = [i for i, s in enumerate(ids) if s in ph.index and np.isfinite(icv[i])]
    sid = [ids[i] for i in keep]
    X, icv = X[keep], icv[keep]
    y   = ph.loc[sid, 'age'].to_numpy(float)
    sex = ph.loc[sid, 'sex'].astype(str).to_numpy()
    ped = ph.loc[sid, 'pedigree'].astype(str).to_numpy()
    return X, icv, vm, y, sex, ped


def main():
    import lightgbm as lgb
    from sklearn.model_selection import StratifiedKFold, GroupKFold

    m = load_primary()
    X, icv, vm, y, sex, ped = load_data(m)
    print(f'N = {len(y)}   {X.shape[1]} features')
    print(f'pre-specified depth grid {m.DEPTH_GRID}, diagnostic grid {DEPTH_GRID}')
    print(f'learning rates {m.LR_GRID}, seed {SEED}, {N_FOLDS} outer fold(s)\n')

    # Exactly the split the primary run uses for this seed.
    fam_strat = m.family_strata(y, sex, ped, m.N_OUTER)
    uped = np.array(sorted(set(ped)))
    ustrat = np.array([fam_strat[g] for g in uped])
    gf = {}
    for f, (_, te) in enumerate(StratifiedKFold(
            m.N_OUTER, shuffle=True, random_state=SEED).split(uped, ustrat)):
        for g in uped[te]:
            gf[g] = f
    fold = np.array([gf[g] for g in ped])

    grid = [(lr, md) for lr in m.LR_GRID for md in DEPTH_GRID]
    rows = []

    for f in range(min(N_FOLDS, m.N_OUTER)):
        tr = np.where(fold != f)[0]
        sc = {c: [] for c in grid}
        for a, b in GroupKFold(m.N_INNER).split(tr, groups=ped[tr]):
            itr, iva = tr[a], tr[b]
            Xi, (Xv,) = m.fit_preprocessing(X, icv, vm, itr, [iva])
            for lr, md in grid:
                mdl = lgb.LGBMRegressor(
                    n_estimators=500, learning_rate=lr, max_depth=md,
                    num_leaves=31, subsample=0.8, subsample_freq=1,
                    colsample_bytree=0.6, verbose=-1,
                    random_state=0).fit(Xi, y[itr])
                sc[(lr, md)].append(m.r2_ref(y[iva], mdl.predict(Xv), y[itr].mean()))
        for (lr, md), v in sc.items():
            rows.append(dict(seed=SEED, outer_fold=f, learning_rate=lr,
                             max_depth=md, inner_r2=float(np.mean(v)),
                             in_prespecified_grid=md in m.DEPTH_GRID))
        best = max(grid, key=lambda c: np.mean(sc[c]))
        print(f'  fold {f}: best lr={best[0]} depth={best[1]} '
              f'(inner R2 {np.mean(sc[best]):.4f})', flush=True)

    out = pd.DataFrame(rows)

    # Averaged over the folds examined. Depth is the question, so the learning
    # rate is collapsed by taking the best one at each depth, which is what the
    # inner loop would have done.
    print('\nInner R2 by depth, best learning rate at each depth, '
          f'mean over {min(N_FOLDS, m.N_OUTER)} outer fold(s):\n')
    print(f"  {'depth':>6}{'inner R2':>12}{'best lr':>10}   {'in grid':<8}")
    summary = []
    for md in DEPTH_GRID:
        sub = out[out.max_depth == md].groupby('learning_rate').inner_r2.mean()
        lr_best = sub.idxmax()
        print(f'  {md:>6}{sub.max():>12.4f}{lr_best:>10}   '
              f'{"yes" if md in m.DEPTH_GRID else "NOT TESTED in the primary run":<8}')
        summary.append((md, sub.max()))

    best_depth = max(summary, key=lambda t: t[1])[0]
    print()
    if best_depth == min(m.DEPTH_GRID):
        print(f'  Depth {best_depth} wins, and it is the lowest value in the')
        print('  pre-specified grid. Depths below it were tested here and did')
        print('  not do better, so the boundary is real and can be reported as')
        print('  checked rather than as an open question.')
    elif best_depth < min(m.DEPTH_GRID):
        d3 = dict(summary)[min(m.DEPTH_GRID)]
        print(f'  Depth {best_depth} beats depth {min(m.DEPTH_GRID)} by '
              f'{max(t[1] for t in summary) - d3:.4f} in inner R2.')
        print('  The pre-specified grid was too narrow at the lower end. Report')
        print('  this as a limitation. Do NOT widen the grid and refit: the')
        print('  reported gaps come from the registered specification.')
    else:
        print(f'  Depth {best_depth} wins, inside the pre-specified grid.')
        print('  Nothing to report beyond the selection counts already logged.')

    os.makedirs(RESD, exist_ok=True)
    p = os.path.join(RESD, 'depth_diagnostic.csv')
    out.to_csv(p, index=False)
    print(f'\nwrote {p}')
    print('Diagnostic only. No primary output was read or written.')


if __name__ == '__main__':
    main()
