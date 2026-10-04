"""
STEP 7, Sensitivity quantities that steps 4 to 6 do not produce.

Adds three things and changes nothing in the primary analysis.

  1. RESIDUALISED BAG
     An alternative age-bias correction: the uncorrected gap regressed on chronological
     age, residual kept. The primary correction instead rescales the prediction with
     cross-fitted coefficients estimated inside each training fold. The two put BAG on
     different scales, so the comparison is about the direction and relative size of the
     downstream coefficients, not about matching numbers.

     Two properties worth keeping in mind when reading the output:
       - The residualisation is fitted on the whole sample, since that is what the
         procedure is: one regression of the gap on age. It is not cross-fitted.
       - Age is already a covariate in H1 and H3, so residualising on age and then
         adjusting for age is close to redundant there. The comparison carries more
         information for H2, where BAG is the outcome.

  2. PCA VERSUS NO-PCA AGREEMENT
     Correlation between the BAG estimates from the reduced and unreduced functional
     models. Neither step 4 nor step 6 can compute it, because each sees only one run.

  3. ACCURACY UNDER EACH CORRECTION
     Pooled R-squared, MAE, RMSE, mean error, slope against chronological age and the
     gap-age correlation, for the uncorrected, cross-fitted and residualised versions
     of every run.

     The pre-specified accuracy of the study is the UNCORRECTED one, which steps 4 and
     5 already write to their metrics files. What is added here is the same quantities
     after correction, because figures and text that show them need a logged source
     rather than a number recomputed inside a plotting script.

     Reading them: the correction divides by a slope below one, which stretches the
     predictions away from the sample mean. That restores the slope against age to one
     and removes the age dependence of the gap, and it necessarily moves predictions
     further from the truth, so R-squared falls and MAE rises. A lower corrected
     R-squared is the expected cost of an unbiased gap, not a worse model.

  4. STABILITY SUMMARIES
     Across-partition SD of each gap and the number of repetitions in which the age-bias
     correction was valid. Reported only; no participant is excluded on these grounds.

INPUTS   bag_functional_pca.csv, bag_functional_nopca.csv, bag_structural.csv
         (whichever exist; missing ones are skipped)
OUTPUT   ../../data/clean/bag_sensitivity.csv
         ../../results/bag_accuracy_by_correction.csv
             participant_id, age, and per run:
             BAG_<name>_corrected, BAG_<name>_raw, BAG_<name>_residualised

The residualised columns are inputs to the R sensitivity models; producing them is not
the sensitivity analysis.
"""
# ============================================================================
# Step 7 of 8. Sensitivity quantities only; changes nothing in the primary analysis.
# Adds an alternative age-bias correction (residualised rather than rescaled) and
# the per-participant stability of the gap across the ten repeated partitions.
# Nothing downstream depends on this script.
# ============================================================================
import os, sys, numpy as np, pandas as pd

WORK = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(WORK, '..', '..'))
OUTD = os.path.join(ROOT, 'data', 'clean')
RESD = os.path.join(ROOT, 'results')

# name -> (filename, BAG column written by the producing script)
RUNS = {
    'functional_pca':   ('bag_functional_pca.csv',   'BAG_functional'),
    'functional_nopca': ('bag_functional_nopca.csv', 'BAG_functional'),
    'structural':       ('bag_structural.csv',       'BAG_structural'),
}


def residualise(gap, age):
    """Regress the gap on chronological age and return the residual.

    Fitted on every participant with a finite value. Participants with a missing gap stay
    missing.
    """
    ok = np.isfinite(gap) & np.isfinite(age)
    out = np.full(len(gap), np.nan)
    if ok.sum() < 3:
        return out, np.nan, np.nan
    slope, intercept = np.polyfit(age[ok], gap[ok], 1)
    out[ok] = gap[ok] - (intercept + slope * age[ok])
    return out, float(slope), float(intercept)


def main():
    frames, loaded = {}, []
    for name, (fname, col) in RUNS.items():
        path = os.path.join(WORK, fname)
        if not os.path.exists(path):
            print(f'  {name}: {fname} not found, skipped')
            continue
        d = pd.read_csv(path)
        need = {'participant_id', 'age', 'pred_age_raw', col}
        if not need.issubset(d.columns):
            sys.stderr.write(f'WARNING: {fname} lacks {sorted(need - set(d.columns))}; '
                             f'{name} skipped.\n')
            continue
        frames[name] = d
        loaded.append(name)
        print(f'  {name}: {len(d)} participants')
    if not frames:
        raise SystemExit('no BAG files found, run steps 4 and 5 first')

    # ---- assemble, and residualise each run -------------------------------------
    base = frames[loaded[0]][['participant_id', 'age']].copy()
    tab = base.copy()
    print('\nAlternative age-bias correction (gap residualised on chronological age)')
    print(f'  {"run":<18} {"slope":>8} {"r(corrected, age)":>19} '
          f'{"r(residualised, age)":>21} {"r(corrected, resid.)":>21}')
    for name in loaded:
        d, col = frames[name], RUNS[name][1]
        d = d.set_index('participant_id').reindex(tab['participant_id'])
        age = d['age'].to_numpy(float)
        corrected = d[col].to_numpy(float)
        raw = d['pred_age_raw'].to_numpy(float) - age          # UNCORRECTED gap
        resid, slope, _ = residualise(raw, age)

        tab[f'BAG_{name}_corrected']    = corrected
        tab[f'BAG_{name}_raw']          = raw
        tab[f'BAG_{name}_residualised'] = resid

        def r(u, v):
            m = np.isfinite(u) & np.isfinite(v)
            return np.corrcoef(u[m], v[m])[0, 1] if m.sum() > 2 else np.nan
        print(f'  {name:<18} {slope:>8.3f} {r(corrected, age):>19.3f} '
              f'{r(resid, age):>21.3f} {r(corrected, resid):>21.3f}')
    print('  A residualised gap is uncorrelated with age by construction, so the third\n'
          '  column checks the arithmetic. The last column shows how far the two\n'
          '  corrections rank participants the same way.')

    # ---- PCA versus no-PCA agreement --------------------------------------------
    if 'functional_pca' in frames and 'functional_nopca' in frames:
        a = tab['BAG_functional_pca_corrected'].to_numpy(float)
        b = tab['BAG_functional_nopca_corrected'].to_numpy(float)
        m = np.isfinite(a) & np.isfinite(b)
        print(f'\nPCA versus unreduced functional BAG (n = {int(m.sum())})')
        print(f'  Pearson r        {np.corrcoef(a[m], b[m])[0, 1]:.3f}')
        print(f'  mean difference  {np.mean(a[m] - b[m]):+.3f} years   '
              f'SD of difference {np.std(a[m] - b[m]):.3f}')
    else:
        print('\nPCA versus unreduced functional BAG: needs both bag_functional_pca.csv '
              'and bag_functional_nopca.csv, skipped')

    # ---- stability summaries ------------------------------------------------------
    print('\nStability across repeated partitions (descriptive; no one is excluded)')
    for name in loaded:
        d = frames[name]
        if 'BAG_sd_across_seeds' in d.columns:
            s = d['BAG_sd_across_seeds'].to_numpy(float)
            s = s[np.isfinite(s)]
            if len(s):
                q = np.percentile(s, [50, 90, 95, 99])
                print(f'  {name:<18} SD across partitions  median {q[0]:.3f}  '
                      f'p90 {q[1]:.3f}  p95 {q[2]:.3f}  p99 {q[3]:.3f}  max {s.max():.3f}')
        if 'BAG_n_valid_repetitions' in d.columns:
            v = d['BAG_n_valid_repetitions'].to_numpy(float)
            print(f'  {name:<18} valid corrected repetitions  min {int(np.nanmin(v))}  '
                  f'mean {np.nanmean(v):.2f}  with none {int((v == 0).sum())}')

    # ---- accuracy under each correction -------------------------------------------
    # Pooled over all participants, matching the pooled metrics in steps 4 and 5 so
    # that the uncorrected rows here reproduce them exactly. A mismatch means the
    # sensitivity table and the primary metrics have drifted apart.
    age = tab['age'].to_numpy(float)
    rows = []
    for name in loaded:
        for kind in ('raw', 'corrected', 'residualised'):
            col = f'BAG_{name}_{kind}'
            if col not in tab.columns:
                continue
            g = tab[col].to_numpy(float)
            ok = np.isfinite(g) & np.isfinite(age)
            a, pred = age[ok], age[ok] + g[ok]
            ss_res = float(((a - pred) ** 2).sum())
            ss_tot = float(((a - a.mean()) ** 2).sum())
            slope = float(np.polyfit(a, pred, 1)[0])
            rows.append(dict(
                run=name, correction=kind, n=int(ok.sum()),
                r2_pooled=1 - ss_res / ss_tot,
                mae=float(np.abs(pred - a).mean()),
                rmse=float(np.sqrt(((pred - a) ** 2).mean())),
                mean_error=float((pred - a).mean()),
                slope_vs_age=slope,
                gap_age_corr=float(np.corrcoef(g[ok], a)[0, 1])))

    acc = pd.DataFrame(rows)
    print('\nAccuracy under each correction (pooled over participants)')
    print('  The pre-specified accuracy of the study is the raw row.')
    print(f"  {'run':<18}{'correction':<15}{'R2':>8}{'MAE':>8}{'RMSE':>8}"
          f"{'mean err':>10}{'slope':>8}{'r(gap,age)':>12}")
    for _, r in acc.iterrows():
        print(f"  {r['run']:<18}{r['correction']:<15}{r['r2_pooled']:>8.3f}"
              f"{r['mae']:>8.2f}{r['rmse']:>8.2f}{r['mean_error']:>10.2f}"
              f"{r['slope_vs_age']:>8.2f}{r['gap_age_corr']:>+12.3f}")

    os.makedirs(RESD, exist_ok=True)
    pa = os.path.join(RESD, 'bag_accuracy_by_correction.csv')
    acc.to_csv(pa, index=False)
    print(f'wrote {pa}')

    p = os.path.join(OUTD, 'bag_sensitivity.csv')
    tab.to_csv(p, index=False)
    print(f'\nwrote {p}   rows {len(tab)}   columns {tab.shape[1]}')
    print('The residualised columns still have to be carried into the R models.')


if __name__ == '__main__':
    main()
