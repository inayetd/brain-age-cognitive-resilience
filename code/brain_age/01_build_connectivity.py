"""
STEP 1 of the functional pipeline.

Assembles the 379x379 resting-state CORRELATION matrix for every baseline (V1) session
from the seed-wise release files, and writes meta.json (region order + session IDs).

The release stores Fisher-z values with the diagonal zeroed. Both are left untouched here
and undone in 04_functional_bag.py.

Outputs (in this directory):
    corr_V1.npy   (N, 379, 379) float32
    meta.json     {"rois": [...379], "ids": [...N]}
"""
# ============================================================================
# Step 1 of 8 in the brain-age pipeline. Run first.
# Reads the released seed-wise resting-state files and assembles one 379x379
# connectivity matrix per participant. No modelling, no exclusions, no statistics:
# this script only reshapes the release into an array the later steps can use.
# Run before: 04_functional_bag.py
# ============================================================================
import os, sys, json, numpy as np, pandas as pd

# --- path to the released correlation folder ----------------------------------
D = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 '..', '..', 'data', 'rest_fmri_correlation')
# ------------------------------------------------------------------------------
OUT = os.path.dirname(os.path.abspath(__file__))
P   = 'rfMRI_REST_FullCorrelationConnectivity_'
norm = lambda s: s.replace('-', '_')      # filenames use '-', column headers use '_'

files   = sorted(f for f in os.listdir(D) if f.startswith(P) and f.endswith('.csv'))
seedmap = {norm(f[len(P):-4]): f for f in files}

# Canonical region order = the column order of any seed file. Only the 379 parcels that
# are BOTH a column and a seed file are used; the release also ships network-level seed
# files, which are aggregates and are excluded automatically by this matching.
cols = pd.read_csv(os.path.join(D, files[0]), nrows=0).columns.tolist()[1:]
rois = [norm(c) for c in cols]
missing = [r for r in rois if r not in seedmap]
if missing:
    raise SystemExit(f'{len(missing)} regions have no seed file, e.g. {missing[:3]}')
if len(rois) != 379:
    raise SystemExit(f'expected 379 regions, found {len(rois)}')

sess = pd.read_csv(os.path.join(D, seedmap[rois[0]]), index_col=0, usecols=[0]).index.tolist()
keep = [i for i, s in enumerate(sess) if str(s).endswith('_V1')]
ids  = [sess[i] for i in keep]
N    = len(ids)

meta_path = os.path.join(OUT, 'meta.json')
if not os.path.exists(meta_path):
    json.dump({'rois': rois, 'ids': ids}, open(meta_path, 'w'))
    print(f'wrote meta.json  ({len(rois)} regions, {N} V1 sessions)')
else:
    m = json.load(open(meta_path))
    if m['rois'] != rois or m['ids'] != ids:
        raise SystemExit('existing meta.json disagrees with these files, delete it and rerun')

p  = os.path.join(OUT, 'corr_V1.npy')
mm = np.lib.format.open_memmap(p, mode=('r+' if os.path.exists(p) else 'w+'),
                               dtype=np.float32, shape=(N, 379, 379))
lo = int(sys.argv[1]) if len(sys.argv) > 2 else 0
hi = min(int(sys.argv[2]), 379) if len(sys.argv) > 2 else 379
for i in range(lo, hi):
    d = pd.read_csv(os.path.join(D, seedmap[rois[i]]), index_col=0)
    d.columns = [norm(c) for c in d.columns]
    mm[:, i, :] = d.iloc[keep].reindex(rois, axis=1).values.astype(np.float32)
    if (i - lo) % 50 == 0:
        print(f'  seed {i}/{hi}', flush=True)
mm.flush()

if hi >= 379:
    # Every matrix is checked, not just the first. The seed files are written
    # independently per region, so an asymmetry means a region was misaligned.
    worst, worst_i = 0.0, -1
    for k in range(N):
        A = np.asarray(mm[k], dtype=np.float64)
        a = float(np.abs(A - A.T).max())
        if a > worst:
            worst, worst_i = a, k
    json.dump({'max_asymmetry': worst, 'worst_session': ids[worst_i], 'n_checked': N},
              open(os.path.join(OUT, 'symmetry_check.json'), 'w'), indent=2)
    print(f'done. N={N}   symmetry check over ALL matrices: max|A-A.T| = {worst:.3g} '
          f'(worst: {ids[worst_i]})')
    if worst > 1e-4:
        sys.stderr.write(f'WARNING: max asymmetry {worst:.3g} exceeds 1e-4. Inspect '
                         f'{ids[worst_i]} before using corr_V1.npy.\n')
else:
    print(f'seeds {lo}-{hi} written; rerun for the remaining seeds before using corr_V1.npy')
