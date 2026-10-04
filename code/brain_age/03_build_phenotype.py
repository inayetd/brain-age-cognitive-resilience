"""
STEP 3, Build the single clean analysis phenotype file.

Merges the three verified sources into one table and writes it to
    data/clean/analysis_phenotype.csv

Every later script reads only that file. Nothing downstream touches the raw release
files again.

SOURCES AND WHY THESE ONES
  data/cog_clean_TMTlog.csv          cognitive indicators AND psychosocial indicators.
                                     Comma-delimited, undamaged. Contains
                                     tlbx_loneliness_raw, which the CFA3 model needs.
  data/buffer_factors_V1_final_.csv  age, sex, education, income. Semicolon-delimited.
                                     Only these columns are taken from it, its PSQI and
                                     IPAQ component columns were damaged by a spreadsheet
                                     round-trip (times and dates in numeric fields).
  AABC_Release2_Non-imaging_Data-XL.csv   pedid, site, and the pre-scored PSQI global and
                                     IPAQ category, taken from the release rather than
                                     from the damaged working copy.

FILES DELIBERATELY NOT USED
  AABC_V1_full_merged.csv            loses participants: 'age_open' is the string
                                     "90 or older" for the oldest, which a naive numeric
                                     conversion turns into NaN.
  AABC_Release2_V1_NIH-Toolbox_rawscores.csv
                                     contains tlbx_sec (self-efficacy) where CFA3
                                     requires loneliness, the wrong seven scales.
  V1_AllScores.csv                   spreadsheet-damaged columns.
"""
# ============================================================================
# Step 3 of 8. The single point where participants are selected.
# Every merge, recode and exclusion happens here and nowhere else, so that the
# measurement models, the brain-age models and the hypothesis tests all describe
# the same people. Downstream scripts never reopen a raw release file.
# Run before: everything else.
# ============================================================================
import os, sys, numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ROOT, 'data')
OUTD = os.path.join(DATA, 'clean')
os.makedirs(OUTD, exist_ok=True)

COG    = os.path.join(DATA, 'cog_clean_TMTlog.csv')
BUF    = os.path.join(DATA, 'buffer_factors_V1_final_.csv')
NONIMG = os.path.join(ROOT, 'AABC_Release2_Non-imaging_Data-XL.csv')

COGNITIVE = {                       # seven indicators of latent cognitive performance
    'flanker':      'tlbx_flanker_uncorr_ss',
    'card_sort':    'tlbx_dccs_uncorr_ss',
    'list_sorting': 'tlbx_listsort_uncorr_ss',
    'pattern_comp': 'tlbx_patterncomp_uncorr_ss',
    'picture_seq':  'tlbx_pictseq_uncorr_ss',
    'ravlt':        'ravlt_immediate_recall',
    # Trail Making B is a completion time, so slower = worse. tmtb_rev is -log(seconds),
    # which orients it with the other six indicators, all scored higher = better. Using
    # tmtb_log would leave one negatively loading indicator in the bifactor model.
    'trails_b_rev': 'tmtb_rev',
}
WELLBEING = {                       # psychological well-being factor
    'positive_affect':   'tlbx_posaffect_raw',
    'life_satisfaction': 'tlbx_lifesatisf_raw',
    'meaning_purpose':   'tlbx_meaning_raw',
}
CONNECTED = {                       # social connectedness factor
    'emotional_support':    'tlbx_emotsupport_raw',
    'instrumental_support': 'tlbx_instsupport_raw',
    'friendship':           'tlbx_friendship_raw',
    'loneliness':           'tlbx_loneliness_raw',      # reverse-scored below
}
STRESS = {'perceived_stress': 'tlbx_stress_raw'}


def main():
    for p in (COG, BUF):
        if not os.path.exists(p):
            raise SystemExit(f'missing required file: {p}')

    # Start from demographics, not cognition: the brain-age models need imaging and age
    # only, so a participant with imaging but no cognitive row must still enter them.
    buf0 = pd.read_csv(BUF, sep=';', low_memory=False).rename(
        columns={'id': 'participant_id', 'age_open': 'age'})
    buf0['participant_id'] = buf0['participant_id'].astype(str)
    out = pd.DataFrame({'participant_id': buf0['participant_id'].unique()})

    cog = pd.read_csv(COG, low_memory=False).rename(columns={'id': 'participant_id'})
    cog['participant_id'] = cog['participant_id'].astype(str)
    missing = [v for v in {**COGNITIVE, **WELLBEING, **CONNECTED, **STRESS}.values()
               if v not in cog.columns]
    if missing:
        raise SystemExit(f'columns absent from cog_clean_TMTlog.csv: {missing}')
    ren = {old: new for new, old in {**COGNITIVE, **WELLBEING, **CONNECTED, **STRESS}.items()}
    cg = cog[['participant_id'] + list(ren)].rename(columns=ren)
    for c in ren.values():
        cg[c] = pd.to_numeric(cg[c], errors='coerce')
    out = out.merge(cg, on='participant_id', how='left')

    # Loneliness is reverse-scored so that all four social indicators point towards
    # greater connectedness. Reversal is around the observed maximum.
    out['loneliness_rev'] = out['loneliness'].max() - out['loneliness']

    out['n_cognitive']  = out[list(COGNITIVE)].notna().sum(axis=1)   # 0 if absent from cog file
    out['sem_eligible'] = out['n_cognitive'] >= 5

    # --- demographics -----------------------------------------------------------
    buf = buf0
    # The oldest participants are stored as the string "90 or older". The flag comes from
    # that string before conversion; testing age == 90 afterwards would also catch
    # participants whose released age is exactly 90.
    age_str = buf['age'].astype(str).str.strip()
    topcoded = age_str == '90 or older'
    n_topcode = int(topcoded.sum())
    buf = buf.assign(
        age_topcoded=topcoded,
        age=pd.to_numeric(age_str.str.replace('90 or older', '90', regex=False),
                          errors='coerce'))
    keep = ['participant_id', 'age', 'age_topcoded', 'sex', 'education', 'dm15', 'income']
    out = out.merge(buf[[c for c in keep if c in buf.columns]], on='participant_id', how='left')
    out = out.rename(columns={'dm15': 'education_years'})
    out['income_log'] = np.log(out['income'].where(out['income'] > 0))

    # --- design variables and release-scored lifestyle measures ------------------
    if os.path.exists(NONIMG):
        nz = pd.read_csv(NONIMG, low_memory=False)
        nz.columns = [c.split(' - ')[0].strip().strip('"') for c in nz.columns]
        # Match the dtype used for participant_id everywhere else. A numeric index here
        # would make every reindex below return NaN without raising, which would silently
        # drop pedigree, site, PSQI, IPAQ and the neurological exclusions.
        nz['id'] = nz['id'].astype(str)
        nz = nz[nz['event'] == 'V1'].set_index('id')
        if not nz.index.isin(out.participant_id).any():
            raise SystemExit('no participant in the release file matches the phenotype '
                             f'ids. Release: {nz.index[0]!r}, phenotype: '
                             f'{out.participant_id.iloc[0]!r}')
        take = {'pedigree': 'pedid', 'site': 'site', 'scanner': 'scanner',
                'psqi_global': 'psqi_global', 'ipaq_category': 'ipaq_category',
                'moca': 'moca_sum', 'apoe': 'APOE_Genotype',
                'rsfmri_pct_complete': 'rfMRI_PctCompl', 'rsfmri_runs': 'rfMRI_Count',
                'mr_qc_code': 'MR_QC_Issue_Codes'}
        numeric = {'psqi_global', 'ipaq_category', 'moca',
                   'rsfmri_pct_complete', 'rsfmri_runs'}
        for new, old in take.items():
            if old in nz.columns:
                v = nz[old].reindex(out.participant_id).to_numpy()
                out[new] = pd.to_numeric(v, errors='coerce') if new in numeric else v
        # PSQI is scored so that HIGHER = worse sleep; reverse so higher = better.
        if 'psqi_global' in out.columns:
            out['sleep_quality'] = out['psqi_global'].max() - out['psqi_global']
        # Participants without a pedigree ID each form their own group; left as NaN they
        # would collapse into one false family.
        miss = out['pedigree'].isna()
        out.loc[miss, 'pedigree'] = ['SOLO_%d' % i for i in range(int(miss.sum()))]
        out['pedigree'] = out['pedigree'].astype(str)
        print(f'  pedigree: {(~miss).sum()} with pedid, {miss.sum()} singletons')
        # Neurological exclusions from the head-injury inventory items released for the
        # full cohort: seizures, stroke or brain haemorrhage, brain tumour. Dementia, MCI,
        # Parkinson's and MS are released for a subsample only and are not used.
        NEURO = {'seizure': 'batlq26', 'stroke_haemorrhage': 'batlq35', 'brain_tumour': 'batlq37'}
        flags = []
        for lbl, col in NEURO.items():
            if col in nz.columns:
                v = pd.to_numeric(nz[col].reindex(out.participant_id), errors='coerce').to_numpy()
                out['neuro_' + lbl] = (v == 1)
                flags.append('neuro_' + lbl)
                print(f'  {lbl:<20} ({col}): {int((v == 1).sum())} flagged')
            else:
                sys.stderr.write(f'WARNING: {col} ({lbl}) not in the release file, '
                                 f'this exclusion cannot be applied.\n')
        out['neuro_exclude'] = out[flags].any(axis=1) if flags else False
        print(f'  excluded on neurological history: {int(out.neuro_exclude.sum())}')

    else:
        # Stop rather than warn. Without this file there are no pedigree groups, so folds
        # cannot be family-grouped and accuracy would be inflated by relatives spanning
        # training and test; and the neurological exclusions cannot be applied, so the
        # analytic sample would be the wrong one. Both failures are invisible in the
        # output table itself, which is why a warning is not enough.
        raise SystemExit(
            f'required file not found: {NONIMG}\n'
            '  It supplies pedigree (family grouping), the neurological-exclusion items, '
            'site, and the pre-scored PSQI and IPAQ.\n'
            '  Place it at that exact path and rerun. Do not proceed without it.')

    if 'neuro_exclude' not in out.columns:
        out['neuro_exclude'] = False

    out = out[(out.age >= 36) & out.age.notna() & ~out.neuro_exclude].reset_index(drop=True)

    p = os.path.join(OUTD, 'analysis_phenotype.csv')
    out.to_csv(p, index=False)

    print(f'\nwrote {p}')
    print(f'  rows {len(out)}   age {out.age.mean():.1f} +/- {out.age.std():.1f} '
          f'[{out.age.min():.0f}-{out.age.max():.0f}]   top-coded {n_topcode}')
    print(f'  SEM-eligible (>=5/7 cognitive): {int(out.sem_eligible.sum())}')
    print('\n  coverage:')
    for c in (list(COGNITIVE) + list(WELLBEING) + list(CONNECTED) + list(STRESS)
              + ['income', 'education', 'sleep_quality', 'ipaq_category']):
        if c in out.columns:
            print(f'    {c:<24} n={int(out[c].notna().sum())}')


if __name__ == '__main__':
    main()
