"""
STEP 2 of the structural pipeline.

Assembles the structural feature matrix for every baseline (V1) session from the released
morphometry tables. No modelling here: the brain-age model is 05_structural_bag.py.

Outputs (in this directory):
    Xstruct.npy       (N, n_features) float64, in the row order of meta.json
    icv.npy           (N,) estimated total intracranial volume
    volmask.npy       (n_features,) bool, True for size-dependent features
    struct_cols.json  feature names
    struct_meta.json  the participant IDs actually matched
"""
# ============================================================================
# Step 2 of 8. Independent of step 1; either may run first.
# Assembles the 1,500 structural features (cortical thickness, surface area, cortical
# volume and a T1w/T2w myelin proxy over 360 Glasser regions, plus 60 subcortical
# volumes) and the intracranial volume used later to remove head-size effects.
# Run before: 05_structural_bag.py
# ============================================================================
import re, numpy as np, pandas as pd, json
import os as _os
HERE=_os.path.dirname(_os.path.abspath(__file__))
D=_os.path.join(HERE,'..','..','data','structural','AABC_Release2_StructuralIDPs')
OUT=HERE
# Session labels are not formatted identically across the released tables. The cortical
# tables use <ID>_V1, while asegstats appends a modality suffix, <ID>_V1_MR. A pattern
# anchored to the end of the string therefore matches the cortical tables and silently
# fails on asegstats, which leaves every subcortical feature and the ICV column missing
# without raising anything. Match the visit tag wherever it occurs instead, and drop it
# together with anything after it.
VISIT = r'_V\d+.*$'
ids=[re.sub(VISIT,'',str(s)) for s in json.load(open(f'{OUT}/meta.json'))['ids']]

def load(f,idcol):
    d=pd.read_csv(f'{D}/{f}.csv',low_memory=False)
    s=d[idcol].astype(str)
    pid=s.str.replace(VISIT,'',regex=True)
    vis=s.str.extract(r'_(V\d+)')[0]
    keep=(vis=='V1').to_numpy()
    d=d.loc[keep].drop(columns=[idcol]).set_index(pid[keep])
    out=d.reindex(ids)
    # Hard stop rather than a silent block of NaN. A format mismatch here is invisible
    # downstream: the features simply arrive empty.
    matched=int(out.notna().any(axis=1).sum())
    if matched==0:
        raise SystemExit(f'{f}.csv: no session matched meta.json. First label seen: '
                         f'{s.iloc[0]!r}; expected something like {ids[0]!r}_V1')
    print(f'  {f:<30} matched {matched}/{len(ids)}')
    return out

blocks={}
for f,tag in [('Cortical_Areal_Thicknesses','TH'),('Cortical_Areal_Surface_Areas','SA'),
              ('Cortical_Areal_Volumes','VOL'),('Cortical_Areal_Myelin','MY')]:
    b=load(f,'x___'); b.columns=[f'{tag}_{c}' for c in b.columns]; blocks[tag]=b
aseg=load('asegstats','Session')
icv=aseg['EstimatedTotalIntraCranialVol'].to_numpy(float)
drop=['EstimatedTotalIntraCranialVol','FS_L_SurfaceHoles','FS_R_SurfaceHoles','SurfaceHoles',
      'BrainSegVol-to-eTIV','MaskVol-to-eTIV']
aseg=aseg.drop(columns=[c for c in drop if c in aseg.columns])
aseg.columns=[f'SUB_{c}' for c in aseg.columns]
blocks['SUB']=aseg

X=pd.concat(blocks.values(),axis=1)
# ICV adjustment applies to volumetric features only; thickness and myelin are not
# head-size scaled.
volmask=np.array([c.startswith(('SA_','VOL_','SUB_')) for c in X.columns])
print('structural features: %d  (TH %d, SA %d, VOL %d, MY %d, SUB %d)'%(
    X.shape[1],*[sum(c.startswith(t+'_') for c in X.columns) for t in ['TH','SA','VOL','MY','SUB']]))
print('ICV-adjusted (volumetric): %d   not adjusted (thickness/myelin): %d'%(volmask.sum(),(~volmask).sum()))
print('NaN cells: %d'%np.isnan(X.to_numpy(float)).sum())
# A whole block arriving empty means an ID-format mismatch in that table, not missing data.
allnan=[c for c in X.columns if X[c].isna().all()]
if allnan:
    raise SystemExit(f'{len(allnan)} features are entirely missing, e.g. {allnan[:5]}. '
                     f'This is an ID or column mismatch, not missing data.')
if not np.isfinite(icv).any():
    raise SystemExit('ICV is missing for every participant; asegstats did not match.')
print('ICV: %d of %d finite'%(int(np.isfinite(icv).sum()),len(icv)))
json.dump({'ids': list(X.index)}, open(f'{OUT}/struct_meta.json','w'))
np.save(f'{OUT}/Xstruct.npy',X.to_numpy(float)); np.save(f'{OUT}/icv.npy',icv)
np.save(f'{OUT}/volmask.npy',volmask)
json.dump(X.columns.tolist(),open(f'{OUT}/struct_cols.json','w'))
