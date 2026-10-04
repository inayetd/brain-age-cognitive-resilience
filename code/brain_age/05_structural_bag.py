"""
Structural brain age gap from the released AABC/HCP-A morphometry.

The cross-validation machinery matches 04_functional_bag.py exactly, so that the H4
comparison between modalities is not confounded by the pipeline. The two arms differ only
in the feature set and in the head-size adjustment:

  no manifold step   the structural features are already a plain matrix, so there is no
                     tangent projection and no Karcher mean
  no PCA             1,500 features against roughly 1,390 participants; PCA is applied to
                     the functional arm only
  ICV adjustment     intracranial volume regressed out of the size-dependent features
                     (surface area, cortical and subcortical volumes) using training-fold
                     coefficients. Thickness and myelin are not head-size scaled.

Preprocessing is refitted wherever a model is evaluated:
  inner fold      imputation means, ICV coefficients and scaler from the inner-training
                  participants, applied to the inner-validation fold
  outer fold      the same chain refitted on the full outer-training sample, applied to
                  the outer-test fold
  cross-fitting   the same chain refitted inside each of the N_CALIB folds used for the
                  age-bias intercept and slope

Hyperparameters are chosen once per outer fold and reused in the cross-fitting step
instead of being reselected there, which avoids a third level of nesting.

INPUTS   Xstruct.npy, icv.npy, volmask.npy, struct_cols.json  from 02_build_structural.py
         meta.json                                            from 01_build_connectivity.py
         ../../data/clean/analysis_phenotype.csv              from 03_build_phenotype.py
OUTPUT   bag_structural.csv, bag_structural_metrics.json
"""
# ============================================================================
# Step 5 of 8. Produces the structural brain-age gap.
# The cross-validation machinery is identical to step 4 by design, so that any
# difference between the two modalities reflects the features rather than the
# pipeline. The only structural-specific step is the head-size adjustment.
# Feeds: H1, H2 and H3 for the structural modality.
# ============================================================================
import json, os, re, sys, warnings, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
import lightgbm as lgb

WORK  = os.path.dirname(os.path.abspath(__file__))
ROOT  = os.path.abspath(os.path.join(WORK, '..', '..'))
PHENO = os.path.join(ROOT, 'data', 'clean', 'analysis_phenotype.csv')

N_OUTER, N_INNER = 10, 5
# Five inner folds, matching the functional arm. The inner loop ranks settings rather
# than estimating performance, and preprocessing is refitted inside every inner fold.
LR_GRID    = (0.03, 0.05, 0.1)
DEPTH_GRID = (3, 5, 7)
N_SEEDS    = 10
N_CALIB    = 5                     # cross-fitting folds for the age-bias coefficients
BETA_MIN   = 0.20                  # slope AT or below this -> fold not corrected


def fit_preprocessing(X, icv, vm, fit_idx, apply_idx_list):
    """Fit imputation, ICV adjustment and standardisation on fit_idx; apply elsewhere.

    Everything estimated from data lives here, so no participant outside fit_idx
    influences it. The caller decides what fit_idx is at each level.

    ICV is regressed out of the size-dependent features only, marked by vm. Thickness and
    myelin are not head-size scaled.
    """
    Xf = X[fit_idx].copy()
    Xa = [X[i].copy() for i in apply_idx_list]

    cm = np.nanmean(Xf, axis=0)                       # imputation: training means only
    Xf = np.where(np.isnan(Xf), cm, Xf)
    Xa = [np.where(np.isnan(x), cm, x) for x in Xa]

    A = np.column_stack([np.ones(len(fit_idx)), icv[fit_idx]])
    B = np.linalg.lstsq(A, Xf[:, vm], rcond=None)[0]
    Xf[:, vm] -= A @ B
    for x, idx in zip(Xa, apply_idx_list):
        x[:, vm] -= np.column_stack([np.ones(len(idx)), icv[idx]]) @ B

    mu, sd = Xf.mean(0), Xf.std(0) + 1e-9
    return (Xf - mu) / sd, [(x - mu) / sd for x in Xa]


def r2_ref(y_true, y_pred, ref_mean):
    """R-squared against an explicit reference mean.

    sklearn's r2_score uses the mean of y_true, i.e. the validation fold's own mean.
    Passing the training-fold mean measures improvement over predicting that mean
    instead. Within one validation fold the denominator is the same for every candidate,
    so the choice does not affect which setting wins.
    """
    ss_res = float(((y_true - y_pred) ** 2).sum())
    ss_tot = float(((y_true - ref_mean) ** 2).sum())
    return 1.0 - ss_res / ss_tot


def family_strata(y, sex, ped, n_outer):
    """One stratification label per family, for the outer split.

    Outer folds are cut over families, so the stratifier must be a family-level quantity.
    Median age and modal sex are used because they do not depend on row order; for
    single-participant families they equal that participant's values.

    The number of age bins drops if any stratum would hold fewer families than there are
    outer folds, since StratifiedKFold cannot balance a class it cannot split. Returns a
    dict keyed by family identifier.
    """
    df = pd.DataFrame(dict(y=y, sex=pd.Series(sex).astype(str), ped=pd.Series(ped).astype(str)))
    g = df.groupby('ped')
    fam_age = g['y'].median()
    fam_sex = g['sex'].agg(lambda s: sorted(s.mode())[0])
    for nb in (10, 5, 4, 2, 1):
        if nb == 1:
            lab = fam_sex
            break
        b = pd.qcut(fam_age, nb, labels=False, duplicates='drop').astype(str)
        lab = fam_sex.str.cat(b)
        if lab.value_counts().min() >= n_outer:
            break
    n_multi = int((g.size() > 1).sum())
    print(f'  families: {len(fam_age)} ({n_multi} with more than one member); '
          f'{lab.nunique()} strata, smallest {lab.value_counts().min()}')
    return lab.to_dict()


def one_run(X, icv, vm, y, ped, fam_strat, seed):
    """One complete nested CV. Returns out-of-fold predicted age and fold labels.

    Preprocessing is refitted in each inner fold, again on the outer-training sample once
    hyperparameters are chosen, and again in each cross-fitting fold. Hyperparameters are
    chosen once per outer fold and reused in the cross-fitting step.
    """
    uped = np.unique(ped)
    ustrat = np.array([fam_strat[g] for g in uped])
    gf = {}
    for f, (_, te) in enumerate(StratifiedKFold(N_OUTER, shuffle=True,
                                                random_state=seed).split(uped, ustrat)):
        for g in uped[te]:
            gf[g] = f
    fold = np.array([gf[g] for g in ped])
    grid = [(lr, md) for lr in LR_GRID for md in DEPTH_GRID]

    pred, log = np.zeros(len(y)), []
    for f in range(N_OUTER):
        tr, te = np.where(fold != f)[0], np.where(fold == f)[0]

        # ---- inner loop: preprocessing refitted inside every inner-training fold ----
        sc = {c: [] for c in grid}
        for a, b in GroupKFold(N_INNER).split(tr, groups=ped[tr]):
            itr, iva = tr[a], tr[b]
            Xi, (Xv,) = fit_preprocessing(X, icv, vm, itr, [iva])
            for lr, md in grid:
                m = lgb.LGBMRegressor(n_estimators=500, learning_rate=lr, max_depth=md,
                                      num_leaves=31, subsample=0.8, subsample_freq=1,
                                      colsample_bytree=0.6, verbose=-1,
                                      random_state=0).fit(Xi, y[itr])
                sc[(lr, md)].append(r2_ref(y[iva], m.predict(Xv), y[itr].mean()))
        lr, md = max(grid, key=lambda c: np.mean(sc[c]))
        bs = float(np.mean(sc[(lr, md)]))
        best = dict(learning_rate=lr, max_depth=md)
        mk = lambda: lgb.LGBMRegressor(n_estimators=500, num_leaves=31, subsample=0.8,
                                       subsample_freq=1, colsample_bytree=0.6,
                                       verbose=-1, random_state=0, **best)

        # ---- outer fold: refit the whole pipeline on the full outer-training sample ----
        Xtr, (Xte,) = fit_preprocessing(X, icv, vm, tr, [te])
        pred[te] = mk().fit(Xtr, y[tr]).predict(Xte)

        # ---- cross-fitted training predictions for the age-bias coefficients ----
        ptr = np.zeros(len(tr))
        for a, b in GroupKFold(N_CALIB).split(tr, groups=ped[tr]):
            ctr, cva = tr[a], tr[b]
            Xc, (Xv,) = fit_preprocessing(X, icv, vm, ctr, [cva])
            ptr[b] = mk().fit(Xc, y[ctr]).predict(Xv)
        beta, alpha = np.polyfit(y[tr], ptr, 1)

        # A fold at or below BETA_MIN is left uncorrected: its predictions still count
        # towards accuracy but contribute to no participant's BAG. Strictly greater than,
        # so equality fails the test. Same rule as the functional arm.
        beta_ok = bool(np.isfinite(beta) and beta > BETA_MIN)
        if not beta_ok:
            sys.stderr.write(
                f'WARNING seed {seed} fold {f}: age-bias slope beta={beta:.4f} is at or below '
                f'{BETA_MIN}. This fold will NOT be age-bias corrected.\n')

        # Both R2 conventions are logged:
        #   outer_r2            reference is the test fold's own mean (sklearn default)
        #   outer_r2_trainmean  reference is the training-fold mean, matching the inner
        #                       loop. Folds are age-stratified, so the two stay close.
        ss_res  = float(((y[te] - pred[te]) ** 2).sum())
        o_r2    = r2_score(y[te], pred[te])
        o_r2_tm = 1.0 - ss_res / float(((y[te] - y[tr].mean()) ** 2).sum())
        log.append(dict(seed=seed, fold=f, n_test=int(len(te)), inner_r2=float(bs),
                        outer_r2=float(o_r2),
                        outer_r2_trainmean=float(o_r2_tm),
                        outer_mae=float(mean_absolute_error(y[te], pred[te])),
                        beta=float(beta), alpha=float(alpha), beta_ok=beta_ok, hp=best))
        print(f'  seed {seed} fold {f}: inner R2={bs:.3f} outer R2={o_r2:.3f} '
              f'optimism={bs - o_r2:+.3f} beta={beta:.3f}{"" if beta_ok else " [NOT CORRECTED]"} '
              f'{best}', flush=True)
    return pred, fold, log


def main():
    # Drop the visit tag and anything after it. Some released tables append a
    # modality suffix (<ID>_V1_MR), so the tag is not always terminal.
    ids = [re.sub(r'_V\d+.*$', '', str(s))
           for s in json.load(open(os.path.join(WORK, 'meta.json')))['ids']]
    X   = np.load(os.path.join(WORK, 'Xstruct.npy'))
    icv = np.load(os.path.join(WORK, 'icv.npy'))
    vm  = np.load(os.path.join(WORK, 'volmask.npy'))

    ph  = pd.read_csv(PHENO).set_index('participant_id')
    keep = [i for i, s in enumerate(ids) if s in ph.index and np.isfinite(icv[i])]
    sid  = [ids[i] for i in keep]
    X, icv = X[keep], icv[keep]
    y    = ph.loc[sid, 'age'].to_numpy(float)
    sex  = ph.loc[sid, 'sex'].astype(str).to_numpy()
    ped  = ph.loc[sid, 'pedigree'].astype(str).to_numpy()
    print(f'N = {len(sid)}   {X.shape[1]} features   '
          f'age {y.mean():.1f} +/- {y.std():.1f} [{y.min():.0f}-{y.max():.0f}]')
    print(f'  ICV-adjusted: {int(vm.sum())}   not adjusted (thickness/myelin): {int((~vm).sum())}')

    fam_strat = family_strata(y, sex, ped, N_OUTER)

    # Resume: a checkpoint is reused only if it describes the same participants in the
    # same order and the same configuration.
    cfg = json.dumps(dict(n_outer=N_OUTER, n_inner=N_INNER, n_calib=N_CALIB,
                          n_seeds=N_SEEDS, lr=list(LR_GRID), depth=list(DEPTH_GRID),
                          beta_min=BETA_MIN), sort_keys=True)
    preds, corrs, logs, start_seed = [], [], [], 0
    ck  = os.path.join(WORK, 'bag_structural_checkpoint.npz')
    ckl = os.path.join(WORK, 'bag_structural_checkpoint_logs.json')
    if os.path.exists(ck) and os.path.exists(ckl):
        z = np.load(ck, allow_pickle=True)
        same_cfg = str(z['cfg']) == cfg if 'cfg' in z.files else False
        if list(map(str, z['sid'])) == list(map(str, sid)) and same_cfg:
            preds = [a for a in z['preds']]
            corrs = [a for a in z['corrs']]
            logs  = json.load(open(ckl))
            start_seed = int(z['seeds_done'])
            print(f'  resuming from checkpoint: seeds 0-{start_seed - 1} complete, '
                  f'continuing at seed {start_seed}', flush=True)

    for seed in range(start_seed, N_SEEDS):
        pred, fold, log = one_run(X, icv, vm, y, ped, fam_strat, seed)
        c = np.full(len(y), np.nan)
        for r in log:
            if not r['beta_ok']:
                continue                       # NaN stays, see the averaging below
            te = np.where(fold == r['fold'])[0]
            c[te] = (pred[te] - r['alpha']) / r['beta']
        preds.append(pred); corrs.append(c); logs += log
        # Checkpoint after every seed; both output files are written only at the end.
        np.savez_compressed(os.path.join(WORK, 'bag_structural_checkpoint.npz'),
                            preds=np.array(preds), corrs=np.array(corrs),
                            y=y, sid=np.array(sid), seeds_done=len(preds), cfg=cfg)
        json.dump(logs, open(os.path.join(WORK, 'bag_structural_checkpoint_logs.json'), 'w'))
        print(f'  [checkpoint written: {len(preds)} of {N_SEEDS} seeds complete]', flush=True)

    pred = np.mean(preds, axis=0)
    C = np.array(corrs)                       # (n_seeds, n_participants), NaN where skipped
    n_valid = np.sum(np.isfinite(C), axis=0)
    # nanmean/nanstd warn on participants with no correctable repetition; that case is
    # handled explicitly below.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        corr = np.nanmean(C, axis=0)
        bag = corr - y
        bag_sd = np.nanstd(C - y, axis=0)
    # One correctable repetition gives nanstd == 0, which would read as perfect
    # stability rather than as undefined.
    bag_sd = np.where(n_valid >= 2, bag_sd, np.nan)
    n_dropped = int((n_valid == 0).sum())
    if n_dropped:
        sys.stderr.write(f'WARNING: {n_dropped} participants have no correctable repetition '
                         f'and therefore no corrected BAG.\n')
    n_skipped_folds = int(sum(not r['beta_ok'] for r in logs))
    if n_skipped_folds:
        sys.stderr.write(f'WARNING: {n_skipped_folds} of {len(logs)} outer folds were not '
                         f'age-bias corrected (beta <= {BETA_MIN}).\n')

    ir  = np.array([r['inner_r2'] for r in logs])
    orr = np.array([r['outer_r2'] for r in logs])
    # Selection optimism compares inner with outer R2, so both sides need the same
    # reference mean. The inner loop scores against the training mean.
    orr_tm = np.array([r['outer_r2_trainmean'] for r in logs])
    om  = np.array([r['outer_mae'] for r in logs])
    # Two dispersions, kept apart. Pooling all N_SEEDS * N_OUTER results into one SD
    # would mix fold-to-fold variation with partition-to-partition variation.
    by_seed = [[r for r in logs if r['seed'] == s] for s in range(N_SEEDS)]
    r2_fold_sd  = [float(np.std([r['outer_r2'] for r in g])) for g in by_seed]
    r2_seed_mn  = [float(np.mean([r['outer_r2'] for r in g])) for g in by_seed]
    mae_fold_sd = [float(np.std([r['outer_mae'] for r in g])) for g in by_seed]
    mae_seed_mn = [float(np.mean([r['outer_mae'] for r in g])) for g in by_seed]

    # Distribution of selected hyperparameters across outer folds.
    hp_counts = {}
    for r in logs:
        k = f"lr{r['hp']['learning_rate']}_d{r['hp']['max_depth']}"
        hp_counts[k] = hp_counts.get(k, 0) + 1
    hp_counts = dict(sorted(hp_counts.items(), key=lambda kv: -kv[1]))

    # Error by age band; a single MAE hides the loss of accuracy at the extremes of
    # the training range.
    bands, edges = {}, [36, 50, 60, 70, 80, 91]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (y >= lo) & (y < hi)
        if m.sum():
            bands[f'{lo}-{hi - 1}'] = dict(n=int(m.sum()),
                                           mae=float(mean_absolute_error(y[m], pred[m])),
                                           mean_signed_error=float((pred[m] - y[m]).mean()),
                                           bag_mean=float(bag[m].mean()))

    met = dict(n=int(len(y)), n_features=int(X.shape[1]), n_seeds=N_SEEDS,
               n_outer=N_OUTER, n_inner=N_INNER,
               inner_r2_mean=float(ir.mean()),
               outer_r2_mean=float(orr.mean()),
               outer_r2_trainmean_mean=float(np.mean([r['outer_r2_trainmean'] for r in logs])),
               outer_r2_sd_across_folds_within_seed=float(np.mean(r2_fold_sd)),
               outer_r2_sd_across_seeds=float(np.std(r2_seed_mn)),
               outer_mae_mean=float(om.mean()),
               outer_mae_sd_across_folds_within_seed=float(np.mean(mae_fold_sd)),
               outer_mae_sd_across_seeds=float(np.std(mae_seed_mn)),
               outer_r2_sd_pooled_do_not_report_as_fold_sd=float(orr.std()),
               hyperparameters_selected=hp_counts,
               error_by_age_band=bands,
               optimism_gap=float(ir.mean() - orr_tm.mean()),
               pooled_r=float(np.corrcoef(pred, y)[0, 1]),
               pooled_r2=float(r2_score(y, pred)),
               pooled_mae=float(mean_absolute_error(y, pred)),
               pooled_rmse=float(np.sqrt(((pred - y) ** 2).mean())),
               mean_error=float((pred - y).mean()),
               bag_age_corr_raw=float(np.corrcoef(pred - y, y)[0, 1]),
               bag_age_corr_corrected=float(np.corrcoef(bag[n_valid > 0], y[n_valid > 0])[0, 1]),
               bag_sd_within_participant_mean=float(np.nanmean(bag_sd)),
               beta_min_threshold=BETA_MIN,
               n_folds_not_corrected=n_skipped_folds,
               n_folds_total=len(logs),
               n_participants_without_corrected_bag=n_dropped,
               valid_repetitions_min=int(n_valid.min()),
               valid_repetitions_mean=float(n_valid.mean()),
               folds=logs)
    print(json.dumps({k: v for k, v in met.items() if k != 'folds'}, indent=2))

    pd.DataFrame(dict(participant_id=sid, age=y, pred_age_raw=pred, pred_age_corr=corr,
                      BAG_structural=bag, BAG_sd_across_seeds=bag_sd,
                      BAG_n_valid_repetitions=n_valid)
                 ).to_csv(os.path.join(WORK, 'bag_structural.csv'), index=False)
    json.dump(met, open(os.path.join(WORK, 'bag_structural_metrics.json'), 'w'), indent=2)


if __name__ == '__main__':
    main()
