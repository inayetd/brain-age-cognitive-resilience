"""
Functional brain age gap from AABC/HCP-A resting-state correlation connectivity.

Pipeline
    correlation matrices -> tanh (undo Fisher z) -> unit diagonal -> positive-definiteness
    check with shrinkage fallback -> tangent space at the training-fold Karcher mean ->
    upper off-diagonal triangle (71,631 features) -> standardise -> PCA (N_PC components)
    -> LightGBM in nested CV -> age-bias correction -> BAG.

Preprocessing is refitted wherever a model is evaluated:
    inner fold      tangent reference, scaler and PCA basis from the inner-training
                    participants, applied to the inner-validation fold
    outer fold      the same chain refitted on the full outer-training sample, applied to
                    the outer-test fold
    cross-fitting   the same chain refitted inside each cross-fitting fold used for the
                    age-bias coefficients

Hyperparameters are chosen once per outer fold and reused in the cross-fitting step
instead of being reselected there, which avoids a third level of nesting.

N_PC is fixed, not tuned; only learning rate and tree depth are searched. USE_PCA = False
runs the unreduced specification on all 71,631 features, everything else unchanged.

INPUTS
    corr_V1.npy   (N, 379, 379) float32       from 01_build_connectivity.py
    meta.json     {"rois": [...379], "ids": [...N]}
    ../../data/clean/analysis_phenotype.csv   from 03_build_phenotype.py

OUTPUTS   TAG is 'pca' or 'nopca'
    bag_functional_TAG.csv   participant_id, age, pred_age_raw, pred_age_corr,
                             BAG_functional, BAG_sd_across_seeds,
                             BAG_n_valid_repetitions, sem_eligible, site.
                             No outer_fold column: with repeated seeds a participant is
                             tested once per seed, not once overall.
    bag_metrics_TAG.json     per-fold and per-seed accuracy, selected hyperparameters,
                             error by age band, Karcher convergence, shrinkage report
"""
# ============================================================================
# Step 4 of 8. Produces the functional brain-age gap.
# This is the script to read if you want to check for leakage. Every quantity
# estimated from data - the tangent-space reference, the scaler, the PCA basis and
# the age-bias coefficients - is fitted inside fit_preprocessing() and inside the
# cross-fitting block, both of which only ever see training participants.
# Feeds: H1, H2 and H3 for the functional modality.
# ============================================================================
import json, os, re, sys, warnings, numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error
import lightgbm as lgb

WORK  = os.path.dirname(os.path.abspath(__file__))
ROOT  = os.path.join(WORK, '..', '..')
CORR  = os.path.join(WORK, 'corr_V1.npy')
META  = os.path.join(WORK, 'meta.json')
PHENO = os.path.join(ROOT, 'data', 'clean', 'analysis_phenotype.csv')

N_OUTER, N_INNER = 10, 5
# The inner loop ranks settings rather than estimating performance, so five folds suffice.
# Preprocessing is refitted inside every inner fold, which dominates runtime.
N_PC       = 100                   # FIXED, not tuned. See the note in the docstring.
LR_GRID    = (0.03, 0.05, 0.1)
DEPTH_GRID = (3, 5, 7)
N_SEEDS    = 10                    # repeated nested CV; predictions averaged
N_CALIB    = 5                     # cross-fitting folds for the age-bias coefficients
KARCHER_IT = 5                     # convergence is logged per fold, see karcher()
# Shrinkage fallback for matrices that fail the positive-definiteness check:
# R* = (1 - lambda) R + lambda I, smallest value that works, per participant.
LAMBDA_GRID = (0.001, 0.005, 0.01, 0.05, 0.10)
BETA_MIN   = 0.20                  # slope AT or below this -> fold not corrected, see below
USE_PCA    = True                  # set False for the pre-specified unreduced sensitivity run
# Output files carry the specification in their name so the two runs do not overwrite
# one another.
TAG = 'pca' if USE_PCA else 'nopca'
# The age floor (36 years) is applied in 03_build_phenotype.py, not here. This script
# takes the phenotype file as given.


# ---------------------------------------------------------------- phenotype
def load_pheno():
    """Read the analysis file written by 03_build_phenotype.py.

    All merging and recoding happens there; nothing downstream touches the raw release
    files.
    """
    if not os.path.exists(PHENO):
        raise SystemExit('run 03_build_phenotype.py first')
    ph = pd.read_csv(PHENO)
    if ph['pedigree'].astype(str).str.startswith('SOLO_').sum() == 0 \
       and ph['pedigree'].nunique() == len(ph):
        sys.stderr.write('WARNING: every participant is their own pedigree group. If the '
                         'non-imaging file was missing when 03 ran, folds are NOT '
                         'pedigree-grouped and accuracy may be inflated.\n')
    return ph.set_index('participant_id')


# ---------------------------------------------------------------- features
def to_correlation(Z):
    """Undo the Fisher z transform and restore the unit diagonal.

    The release stores Fisher-z values with the diagonal zeroed. No eigenvalue repair
    happens here; positive definiteness is checked separately in main().
    """
    A = (Z + Z.T) / 2.0
    R = np.tanh(A)
    np.fill_diagonal(R, 1.0)
    return R


def karcher(S_train, n_iter=KARCHER_IT):
    """Riemannian (Karcher) mean of the training matrices.

    Returns (M^(-1/2), convergence trace). Fixed-point iteration: at the centre the mean
    tangent displacement is zero, so the Frobenius norm of that displacement is the
    convergence criterion. The trace is logged; if the final norm is not small relative
    to the first, KARCHER_IT is too low.
    """
    M = S_train.mean(axis=0).astype(np.float64)
    trace = []
    for _ in range(n_iter):
        w, V = np.linalg.eigh(M); w = np.clip(w, 1e-12, None)
        H, Hi = (V * np.sqrt(w)) @ V.T, (V * w ** -0.5) @ V.T
        acc = np.zeros_like(M)
        for A in S_train:
            C = Hi @ A.astype(np.float64) @ Hi
            e, U = np.linalg.eigh((C + C.T) / 2)
            acc += (U * np.log(np.clip(e, 1e-12, None))) @ U.T
        Sb = acc / len(S_train)
        trace.append(float(np.linalg.norm(Sb)))
        e, U = np.linalg.eigh((Sb + Sb.T) / 2)
        M = H @ ((U * np.exp(e)) @ U.T) @ H; M = (M + M.T) / 2
    w, V = np.linalg.eigh(M)
    return (V * np.clip(w, 1e-12, None) ** -0.5) @ V.T, trace


def tangent(S, Wi, iu):
    out = np.empty((len(S), iu[0].size), dtype=np.float32)
    for i, A in enumerate(S):
        C = Wi @ A.astype(np.float64) @ Wi
        e, U = np.linalg.eigh((C + C.T) / 2)
        out[i] = ((U * np.log(np.clip(e, 1e-12, None))) @ U.T)[iu]
    return out


# ---------------------------------------------------------------- preprocessing block
def fit_preprocessing(S_fit, S_apply_list, iu):
    """Fit the data-dependent preprocessing on S_fit and apply it to S_apply_list.

    The tangent reference, the scaler and the PCA basis are all estimated here and
    nowhere else, so no participant outside S_fit influences them. The caller decides
    what S_fit is at each level.

    Returns (Z_fit, [Z_apply, ...], karcher_trace, explained_variance_ratio). Scores are
    already projected onto the N_PC basis; with USE_PCA False the standardised tangent
    features are returned unreduced and evr is None.
    """
    Wi, ktrace = karcher(S_fit)
    Xf = tangent(S_fit, Wi, iu)
    Xa = [tangent(Sa, Wi, iu) for Sa in S_apply_list]
    mu, sd = Xf.mean(0), Xf.std(0) + 1e-9
    Xf = (Xf - mu) / sd
    Xa = [(x - mu) / sd for x in Xa]
    if not USE_PCA:
        return Xf, Xa, ktrace, None
    _, sv, Vt = np.linalg.svd(Xf, full_matrices=False)
    # Variance share of the retained components. Logged only; N_PC is fixed and is not
    # adjusted on the basis of this value.
    evr = float((sv[:N_PC] ** 2).sum() / (sv ** 2).sum())
    Vt = Vt[:N_PC]
    return Xf @ Vt.T, [x @ Vt.T for x in Xa], ktrace, evr


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


# ---------------------------------------------------------------- one CV run
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


def one_run(S, y, ped, fam_strat, iu, seed):
    """One complete nested CV. Returns out-of-fold predicted age and fold labels.

    Preprocessing is refitted in each inner fold, again on the outer-training sample once
    hyperparameters are chosen, and again in each cross-fitting fold. Hyperparameters are
    the exception: chosen once per outer fold and reused in the cross-fitting step, which
    avoids a third level of nesting.
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
            Ztr_i, (Zva_i,), _, _ = fit_preprocessing(S[itr], [S[iva]], iu)
            for lr, md in grid:
                m = lgb.LGBMRegressor(n_estimators=500, learning_rate=lr, max_depth=md,
                                      num_leaves=31, subsample=0.8, subsample_freq=1,
                                      colsample_bytree=0.6, verbose=-1,
                                      random_state=0).fit(Ztr_i, y[itr])
                sc[(lr, md)].append(r2_ref(y[iva], m.predict(Zva_i), y[itr].mean()))
            del Ztr_i, Zva_i
        lr, md = max(grid, key=lambda c: np.mean(sc[c]))
        bs = float(np.mean(sc[(lr, md)]))
        best = dict(n_pc=(N_PC if USE_PCA else None), learning_rate=lr, max_depth=md)
        hp = dict(learning_rate=lr, max_depth=md)
        mk = lambda: lgb.LGBMRegressor(n_estimators=500, num_leaves=31, subsample=0.8,
                                       subsample_freq=1, colsample_bytree=0.6,
                                       verbose=-1, random_state=0, **hp)

        # ---- outer fold: refit the whole pipeline on the full outer-training sample ----
        Ztr, (Zte,), ktrace, evr = fit_preprocessing(S[tr], [S[te]], iu)
        pred[te] = mk().fit(Ztr, y[tr]).predict(Zte)

        # ---- cross-fitted training predictions for the age-bias coefficients ----
        # Each calibration prediction comes from a pipeline, preprocessing included,
        # fitted without that participant. In-sample predictions would under-correct,
        # being less compressed than out-of-sample ones.
        ptr = np.zeros(len(tr))
        for a, b in GroupKFold(N_CALIB).split(tr, groups=ped[tr]):
            ctr, cva = tr[a], tr[b]
            Zc, (Zv,), _, _ = fit_preprocessing(S[ctr], [S[cva]], iu)
            ptr[b] = mk().fit(Zc, y[ctr]).predict(Zv)
            del Zc, Zv
        beta, alpha = np.polyfit(y[tr], ptr, 1)

        # The correction divides by beta, so a slope near zero turns small prediction
        # differences into very large gaps. A fold at or below BETA_MIN is left
        # uncorrected: its predictions still count towards accuracy but contribute to no
        # participant's BAG. Strictly greater than, so equality fails the test.
        beta_ok = bool(np.isfinite(beta) and beta > BETA_MIN)
        if not beta_ok:
            sys.stderr.write(
                f'WARNING seed {seed} fold {f}: age-bias slope beta={beta:.4f} is at or below '
                f'{BETA_MIN}. This fold will NOT be age-bias corrected. Its predictions still '
                f'contribute to accuracy but not to any participant BAG.\n')

        # Both R2 conventions are logged:
        #   outer_r2            reference is the test fold's own mean (sklearn default)
        #   outer_r2_trainmean  reference is the training-fold mean, matching the inner
        #                       loop. Folds are age-stratified, so the two stay close.
        ss_res  = float(((y[te] - pred[te]) ** 2).sum())
        o_r2    = r2_score(y[te], pred[te])
        o_r2_tm = 1.0 - ss_res / float(((y[te] - y[tr].mean()) ** 2).sum())
        log.append(dict(seed=seed, fold=f, n_test=int(len(te)), inner_r2=float(bs),
                        outer_r2=float(o_r2), outer_mae=float(mean_absolute_error(y[te], pred[te])),
                        outer_r2_trainmean=float(o_r2_tm),
                        beta=float(beta), alpha=float(alpha), beta_ok=beta_ok, hp=best,
                        pca_explained_variance=evr,
                        karcher_trace=ktrace))
        print(f'  seed {seed} fold {f}: inner R2={bs:.3f} outer R2={o_r2:.3f} '
              f'optimism={bs - o_r2:+.3f} beta={beta:.3f}{"" if beta_ok else " [NOT CORRECTED]"} '
              f'karcher {ktrace[0]:.3f}->{ktrace[-1]:.4f} {best}', flush=True)
    return pred, fold, log


# ---------------------------------------------------------------- main
def main():
    meta = json.load(open(META))
    # Drop the visit tag and anything after it. Some released tables append a
    # modality suffix (<ID>_V1_MR), so the tag is not always terminal.
    ids  = [re.sub(r'_V\d+.*$', '', str(s)) for s in meta['ids']]
    raw  = np.load(CORR, mmap_mode='r')

    ph   = load_pheno()
    keep = [i for i, s in enumerate(ids) if s in ph.index]
    sid  = [ids[i] for i in keep]
    y    = ph.loc[sid, 'age'].to_numpy(float)
    sex  = ph.loc[sid, 'sex'].astype(str).to_numpy()
    ped  = ph.loc[sid, 'pedigree'].astype(str).to_numpy()
    print(f'N = {len(sid)}   age {y.mean():.1f} +/- {y.std():.1f}  [{y.min():.0f}-{y.max():.0f}]')

    # The tangent projection takes a matrix logarithm, which needs strictly positive
    # eigenvalues. Shrinkage towards the identity pulls each eigenvalue towards one and
    # preserves their order, so the smallest lambda that clears zero perturbs the matrix
    # least. Repairs and exclusions are recorded in the metrics file.
    S = np.empty((len(keep), 379, 379), dtype=np.float32)
    shrunk, unfixable, ok = {}, [], 0
    for n, i in enumerate(keep):
        R = to_correlation(np.asarray(raw[i], dtype=np.float64))
        if np.linalg.eigvalsh(R).min() > 0:
            ok += 1
            S[n] = R
            continue
        for lam in LAMBDA_GRID:
            Rs = (1.0 - lam) * R + lam * np.eye(379)
            if np.linalg.eigvalsh(Rs).min() > 0:
                shrunk[sid[n]] = lam
                S[n] = Rs
                break
        else:
            unfixable.append(sid[n])
            S[n] = R                      # placeholder; this row is dropped below

    print(f'  positive definite as released: {ok} of {len(keep)}')
    if shrunk:
        by_lam = {l: sum(1 for v in shrunk.values() if v == l) for l in LAMBDA_GRID}
        print(f'  repaired by shrinkage: {len(shrunk)}  ' +
              '  '.join(f'lambda={l}: {c}' for l, c in by_lam.items() if c))
    if unfixable:
        sys.stderr.write(f'WARNING: {len(unfixable)} matrices remain non-positive-definite at '
                         f'lambda={max(LAMBDA_GRID)} and are excluded from the functional '
                         f'analysis: {", ".join(unfixable)}\n')
        drop = {t for t in unfixable}
        m = np.array([t not in drop for t in sid])
        S, y, sex, ped = S[m], y[m], sex[m], ped[m]
        sid = [t for t in sid if t not in drop]
        print(f'  N after exclusions = {len(sid)}')

    pd_report = dict(n_positive_definite_as_released=ok,
                     n_repaired_by_shrinkage=len(shrunk),
                     lambda_grid=list(LAMBDA_GRID),
                     lambda_used=shrunk,
                     n_excluded_not_repairable=len(unfixable),
                     excluded_participants=unfixable)
    iu = np.triu_indices(379, 1)

    fam_strat = family_strata(y, sex, ped, N_OUTER)

    # Repeated nested CV. Averaging over seeds removes the split-dependent part of the
    # participant-level noise; the SD across seeds measures stability.
    #
    # Resume: a checkpoint is reused only if it describes the same participants in the
    # same order AND the same configuration. one_run(seed=k) depends on nothing but k, so
    # a resumed run reproduces an uninterrupted one exactly. The configuration is
    # fingerprinted below because a changed grid or component count would otherwise mix
    # two different analyses in one output file.
    cfg = json.dumps(dict(tag=TAG, n_pc=(N_PC if USE_PCA else None), use_pca=USE_PCA,
                          n_outer=N_OUTER, n_inner=N_INNER, n_calib=N_CALIB,
                          n_seeds=N_SEEDS, lr=list(LR_GRID), depth=list(DEPTH_GRID),
                          karcher_it=KARCHER_IT, beta_min=BETA_MIN,
                          lam=list(LAMBDA_GRID)), sort_keys=True)
    preds, corrs, logs, start_seed = [], [], [], 0
    ck  = os.path.join(WORK, f'bag_functional_{TAG}_checkpoint.npz')
    ckl = os.path.join(WORK, f'bag_functional_{TAG}_checkpoint_logs.json')
    if os.path.exists(ck) and os.path.exists(ckl):
        z = np.load(ck, allow_pickle=True)
        same_cfg = str(z['cfg']) == cfg if 'cfg' in z.files else False
        if list(map(str, z['sid'])) == list(map(str, sid)) and same_cfg:
            preds = [a for a in z['preds']]
            corrs = [a for a in z['corrs']]
            logs  = json.load(open(ckl))
            start_seed = int(z['seeds_done'])
            print(f'  resuming from checkpoint: seeds 0-{start_seed - 1} already complete, '
                  f'continuing at seed {start_seed}', flush=True)
        else:
            sys.stderr.write('WARNING: the checkpoint describes a different participant '
                             'set or a different configuration and was ignored. Starting '
                             'from seed 0.\n')

    for seed in range(start_seed, N_SEEDS):
        pred, fold, log = one_run(S, y, ped, fam_strat, iu, seed)
        # NaN marks "not correctable in this repetition" and propagates to the nanmean
        # and the per-participant count below.
        c = np.full(len(y), np.nan)
        for r in log:
            if not r['beta_ok']:
                continue
            te = np.where(fold == r['fold'])[0]
            c[te] = (pred[te] - r['alpha']) / r['beta']
        preds.append(pred); corrs.append(c); logs += log

        # Checkpoint after every seed. This loop runs for the better part of a day and
        # both output files are written only at the end, so an interruption at seed eight
        # would otherwise discard eight seeds of work.
        np.savez_compressed(os.path.join(WORK, f'bag_functional_{TAG}_checkpoint.npz'),
                            preds=np.array(preds), corrs=np.array(corrs),
                            y=y, sid=np.array(sid), seeds_done=len(preds), cfg=cfg)
        json.dump(logs, open(os.path.join(WORK, f'bag_functional_{TAG}_checkpoint_logs.json'), 'w'))
        print(f'  [checkpoint written: {len(preds)} of {N_SEEDS} seeds complete]', flush=True)
    pred = np.mean(preds, axis=0)
    C = np.array(corrs)                       # (n_seeds, n_participants), NaN where skipped
    n_valid = np.sum(np.isfinite(C), axis=0)  # correctable repetitions per participant
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
    # No correctable repetition means no corrected BAG; the count is reported so these
    # participants can be excluded downstream.
    n_dropped = int((n_valid == 0).sum())
    if n_dropped:
        sys.stderr.write(f'WARNING: {n_dropped} participants have no correctable repetition '
                         f'and therefore no corrected BAG. They are written with a missing '
                         f'BAG and must be excluded from analyses requiring it.\n')
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
    #   within_seed  : SD across the outer folds, averaged over seeds
    #   between_seed : SD of the seed means
    by_seed = [[r for r in logs if r['seed'] == s] for s in range(N_SEEDS)]
    r2_fold_sd = [float(np.std([r['outer_r2'] for r in g])) for g in by_seed]
    r2_seed_mn = [float(np.mean([r['outer_r2'] for r in g])) for g in by_seed]
    mae_fold_sd = [float(np.std([r['outer_mae'] for r in g])) for g in by_seed]
    mae_seed_mn = [float(np.mean([r['outer_mae'] for r in g])) for g in by_seed]

    # Selected hyperparameters across all outer folds. One setting winning nearly every
    # fold means the grid is well centred; an even spread means the criterion is flat.
    hp_counts = {}
    for r in logs:
        k = f"lr{r['hp']['learning_rate']}_d{r['hp']['max_depth']}"
        hp_counts[k] = hp_counts.get(k, 0) + 1
    hp_counts = dict(sorted(hp_counts.items(), key=lambda kv: -kv[1]))

    # Error by age band. Brain-age models are least accurate at the extremes of the
    # training range, and a single MAE hides that.
    bands, edges = {}, [36, 50, 60, 70, 80, 91]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (y >= lo) & (y < hi)
        if m.sum():
            bands[f'{lo}-{hi - 1}'] = dict(
                n=int(m.sum()),
                mae=float(mean_absolute_error(y[m], pred[m])),
                mean_signed_error=float((pred[m] - y[m]).mean()),
                bag_mean=float(bag[m].mean()))

    karcher_first = float(np.mean([r['karcher_trace'][0] for r in logs]))
    karcher_last  = float(np.mean([r['karcher_trace'][-1] for r in logs]))

    met = dict(n=int(len(y)), n_seeds=N_SEEDS, n_outer=N_OUTER, n_inner=N_INNER,
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
               karcher_mean_displacement_first_iter=karcher_first,
               karcher_mean_displacement_last_iter=karcher_last,
               optimism_gap=float(ir.mean() - orr_tm.mean()),
               pooled_r=float(np.corrcoef(pred, y)[0, 1]),
               pooled_r2=float(r2_score(y, pred)),
               pooled_mae=float(mean_absolute_error(y, pred)),
               pooled_rmse=float(np.sqrt(((pred - y) ** 2).mean())),
               mean_error=float((pred - y).mean()),
               bag_age_corr_raw=float(np.corrcoef(pred - y, y)[0, 1]),
               bag_age_corr_corrected=float(np.corrcoef(bag[n_valid > 0], y[n_valid > 0])[0, 1]),
               bag_sd_within_participant_mean=float(np.nanmean(bag_sd)),
               pca_explained_variance_mean=(
                   float(np.mean([r['pca_explained_variance'] for r in logs]))
                   if USE_PCA else None),
               positive_definiteness=pd_report,
               beta_min_threshold=BETA_MIN,
               n_folds_not_corrected=n_skipped_folds,
               n_folds_total=len(logs),
               n_participants_without_corrected_bag=n_dropped,
               valid_repetitions_min=int(n_valid.min()),
               valid_repetitions_mean=float(n_valid.mean()),
               folds=logs)
    print(json.dumps({k: v for k, v in met.items() if k != 'folds'}, indent=2))

    # Looked up defensively: this runs after hours of computation, and a KeyError on a
    # descriptive column would discard all of it.
    tab = dict(participant_id=sid, age=y, pred_age_raw=pred, pred_age_corr=corr,
               BAG_functional=bag, BAG_sd_across_seeds=bag_sd,
               BAG_n_valid_repetitions=n_valid)
    for col in ('sem_eligible', 'site'):
        if col in ph.columns:
            tab[col] = ph.loc[sid, col].to_numpy()
        else:
            sys.stderr.write(f"WARNING: '{col}' not in the phenotype file; omitted.\n")
    pd.DataFrame(tab).to_csv(os.path.join(WORK, f'bag_functional_{TAG}.csv'), index=False)
    json.dump(met, open(os.path.join(WORK, f'bag_metrics_{TAG}.json'), 'w'), indent=2)


if __name__ == '__main__':
    main()
