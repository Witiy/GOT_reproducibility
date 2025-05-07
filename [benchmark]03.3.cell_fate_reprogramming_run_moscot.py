import warnings
import moscot as mt
import moscot.plotting as mpl
from moscot.problems.time import TemporalProblem
import scanpy as sc
import pandas as pd
import numpy as np
print(mt.__version__)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FutureWarning)

adata = sc.read_h5ad("../pygot_data/03_cellfate/adata_repro.h5ad")
embedding_key = 'X_pca'
time_key = 'time_info'
ts = np.sort(np.unique(adata.obs[time_key]))
ts_map = dict(zip(ts, range(len(ts))))
adata.obs['time_label'] = adata.obs[time_key].replace(ts_map).astype(int)
time_key = 'time_label'
adata.obs[time_key] = pd.Categorical(adata.obs[time_key].astype(int))

tp0 = TemporalProblem(adata)
tp0 = tp0.prepare(time_key=time_key, joint_attr=embedding_key)
moscot_res = []

for e in [1e-1, 5e-1, 1e-2, 5e-2, 1e-3, 5e-3, 1e-4]:

    tp0 = tp0.solve(epsilon=e, scale_cost="mean", max_iterations=1e7, device='cpu')
    moscot_pred = []
    
    tp0.pull(source=0, target=2, data='state_info', subset=['Reprogrammed'])
    adata.obs['pull_Reprogrammed'] = adata.obs['pull'].copy()
    tp0.pull(source=0, target=2, data='state_info', subset=['Failed'])
    adata.obs['pull_Failed'] = adata.obs['pull'].copy()
    adata.obs['moscot_fate_bias'] = adata.obs['pull_Reprogrammed'] / (adata.obs['pull_Reprogrammed'] + adata.obs['pull_Failed'])
    moscot_pred.append(adata.obs['moscot_fate_bias'].loc[adata.obs.loc[adata.obs[time_key] == 0].index].copy())
    
    tp0.pull(source=1, target=2, data='state_info', subset=['Reprogrammed'])
    adata.obs['pull_Reprogrammed'] = adata.obs['pull'].copy()
    tp0.pull(source=1, target=2, data='state_info', subset=['Failed'])
    adata.obs['pull_Failed'] = adata.obs['pull'].copy()
    adata.obs['moscot_fate_bias'] = adata.obs['pull_Reprogrammed'] / (adata.obs['pull_Reprogrammed'] + adata.obs['pull_Failed'])
    moscot_pred.append(adata.obs['moscot_fate_bias'].loc[adata.obs.loc[adata.obs[time_key] == 1].index].copy())
    moscot_pred = pd.concat(moscot_pred)
    moscot_pred.loc[np.isnan(moscot_pred)] = 0.5
    pd.DataFrame(moscot_pred).to_csv('../pygot_data/03_cellfate/reprogram/moscot_v4.0_{}.csv'.format(e))
