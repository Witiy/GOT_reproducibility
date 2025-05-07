import warnings
import moscot as mt
import moscot.plotting as mpl
from moscot.problems.time import TemporalProblem
import scanpy as sc
import pandas as pd
print(mt.__version__)
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", FutureWarning)
epsilon=1e-3
adata = sc.read('/storage/xuruihong/repository/pygot_data/03_cellfate/Oren/adata.h5ad')
embedding_key = 'X_pca'
time_key = 'stage_numeric'


adata.obs[time_key] = pd.Categorical(adata.obs[time_key])
tp0 = TemporalProblem(adata)
tp0 = tp0.prepare(time_key=time_key, joint_attr=embedding_key)
moscot_res = []
for e in [1e-1, 5e-1, 1e-2, 5e-2, 1e-3, 5e-3, 1e-4]:
    tp0 = tp0.solve(epsilon=e, scale_cost="mean", max_iterations=1e7)
    tp0.pull(source=2, target=3, data='state_info', subset=['14_high'])
    adata.obs['pull_high'] = adata.obs['pull'].copy()
    tp0.pull(source=2, target=3, data='state_info', subset=['14_low'])
    adata.obs['pull_low'] = adata.obs['pull'].copy()
    adata.obs['moscot_fate_bias'] = adata.obs['pull_low'] / (adata.obs['pull_high'] + adata.obs['pull_low'])
    moscot_res.append(adata.obs.copy())
    moscot_res[-1].to_csv('../pygot_data/03_cellfate/Oren/moscot_v4.0_{}.csv'.format(e))
#tp0 = tp0.solve(epsilon=epsilon, scale_cost="mean", max_iterations=1e7)
#tp0.pull(source=2, target=3, data='state_info', subset=['14_high'])
#adata.obs['pull_high'] = adata.obs['pull'].copy()
#tp0.pull(source=2, target=3, data='state_info', subset=['14_low'])
#adata.obs['pull_low'] = adata.obs['pull'].copy()
#adata.obs['moscot_fate_bias'] = adata.obs['pull_low'] / (adata.obs['pull_high'] + adata.obs['pull_low'])
#adata.obs.to_csv('/storage/xuruihong/repository/pygot_data/03_cellfate/Oren/moscot_{}'.format(epsilon))
