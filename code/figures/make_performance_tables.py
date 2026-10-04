"""
Two tables for the thesis, built from the logged metrics rather than typed.

MAIN TABLE  results/table_performance.csv
    What the age-prediction model achieves, and what the gap it produces looks
    like. Accuracy comes from the UNCORRECTED out-of-sample predictions, which
    is the convention in this literature: the correction uses chronological age,
    so anything computed after it describes the corrected quantity rather than
    the model. Accuracy is the mean over the 100 outer folds, with the two
    sources of variability kept apart.

    Beside it, the two numbers that describe the gap actually used downstream:
    its correlation with chronological age before and after correction, which is
    how the correction is shown to have worked, and its standard deviation,
    which is the scale every standardised coefficient is expressed in.

CALIBRATION DIAGNOSTIC  results/table_calibration_diagnostic.csv
    Accuracy of the CORRECTED predicted age. This is not the accuracy of the
    study and is kept in its own table for that reason. It is reported pooled
    over all participants, because the corrected predictions were never scored
    fold by fold, and the level is stated in the table itself so the two cannot
    be mistaken for one another.

    It is worth showing because it makes the price of the correction explicit:
    dividing by a slope below one stretches the predictions away from the sample
    mean, which restores the slope against age to one and necessarily moves the
    predictions further from the truth.

INPUTS   results/structural/bag_structural_metrics.json
         results/functional/bag_metrics_pca.json
         results/bag_sensitivity.csv
"""
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
RES = os.path.join(ROOT, 'results')

RUNS = {
    'Structural':       os.path.join(RES, 'structural', 'bag_structural_metrics.json'),
    'Functional (PCA)': os.path.join(RES, 'functional', 'bag_metrics_pca.json'),
}
KEY = {'Structural': 'structural', 'Functional (PCA)': 'functional_pca'}


def main():
    sens = pd.read_csv(os.path.join(RES, 'bag_sensitivity.csv'))
    age = sens['age'].to_numpy(float)

    main_rows, diag_rows = [], []
    for label, path in RUNS.items():
        if not os.path.exists(path):
            raise SystemExit(f'not found: {path}')
        m = json.load(open(path))
        k = KEY[label]
        gap = sens[f'BAG_{k}_corrected'].to_numpy(float)
        ok = np.isfinite(gap)

        main_rows.append({
            'Run': label,
            'n': m['n'],
            'R2 (outer fold mean)': round(m['outer_r2_mean'], 3),
            'R2 SD across folds': round(m['outer_r2_sd_across_folds_within_seed'], 3),
            'R2 SD across partitions': round(m['outer_r2_sd_across_seeds'], 4),
            'MAE years (outer fold mean)': round(m['outer_mae_mean'], 2),
            'MAE SD across folds': round(m['outer_mae_sd_across_folds_within_seed'], 2),
            'MAE SD across partitions': round(m['outer_mae_sd_across_seeds'], 3),
            'gap-age r uncorrected': round(m['bag_age_corr_raw'], 3),
            'gap-age r corrected': round(m['bag_age_corr_corrected'], 3),
            'BAG SD years': round(float(gap[ok].std(ddof=1)), 2),
            'BAG within-participant SD': round(m['bag_sd_within_participant_mean'], 2),
        })

        pred = age[ok] + gap[ok]
        a = age[ok]
        diag_rows.append({
            'Run': label,
            'Level': 'pooled over participants',
            'Quantity': 'corrected predicted age',
            'R2': round(1 - float(((a - pred) ** 2).sum())
                        / float(((a - a.mean()) ** 2).sum()), 3),
            'MAE years': round(float(np.abs(pred - a).mean()), 2),
            'RMSE years': round(float(np.sqrt(((pred - a) ** 2).mean())), 2),
            'slope vs chronological age': round(float(np.polyfit(a, pred, 1)[0]), 2),
        })

    main = pd.DataFrame(main_rows)
    diag = pd.DataFrame(diag_rows)

    print('MAIN PERFORMANCE TABLE')
    print('  accuracy from the uncorrected predictions, outer-fold means\n')
    print(main.to_string(index=False))
    print('\n\nCALIBRATION DIAGNOSTIC')
    print('  accuracy of the corrected predicted age; not the accuracy of the study\n')
    print(diag.to_string(index=False))

    main.to_csv(os.path.join(RES, 'table_performance.csv'), index=False)
    diag.to_csv(os.path.join(RES, 'table_calibration_diagnostic.csv'), index=False)
    print(f"\nwrote {os.path.join('results', 'table_performance.csv')}")
    print(f"wrote {os.path.join('results', 'table_calibration_diagnostic.csv')}")


if __name__ == '__main__':
    main()
