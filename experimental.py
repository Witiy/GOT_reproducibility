
import numpy as np
import pandas as pd
import scanpy as sc

def split_time_series(adata, time_key='PseudoTime', k=4):
    max_t = np.max(adata.obs[time_key])
    adata.obs['Day'] = 0
    for i in range(k):
        adata.obs.loc[(adata.obs[time_key] >= i * max_t / k) & (adata.obs[time_key] < (i + 1)* max_t / k), 'Day'] = i
    adata.obs.loc[adata.obs[time_key] == max_t, 'Day'] = k-1
    return adata
def construct_adata(data, obs, ref_network, time_key='PseudoTime', split_k = 5):
    adata = sc.AnnData(data.T, obs=obs)
    adata.obsm['X_origin'] = adata.X
    adata = split_time_series(adata, time_key, k=split_k)
    #tsner = TSNE(perplexity=100)
    #adata.obsm['X_tsne'] = tsner.fit_transform(adata.obsm['X_pca'])
    adata.var['idx'] = range(len(adata.var))
    TF_names = np.unique(ref_network.Gene1)
    TF_idx = adata.var.loc[TF_names]['idx'].to_numpy()
    return adata, TF_names, TF_idx



def experiment(path, split_k = 4):
    data = pd.read_csv(path+'/ExpressionData.csv', index_col=0)
    obs = pd.read_csv(path+'/PseudoTime.csv', index_col=0).loc[data.columns]
    obs['Time'] = obs.index.to_series().apply(lambda x: int(x.split('_')[-1]))
    #cluster = pd.read_csv(path+'/ClusterIds.csv', index_col=0)
    ref_network = pd.read_csv(path+'/refNetwork.csv')

    #data = data.loc[data.index[data.index.to_series().apply(lambda x: int(x[1:])).argsort()]]
    #obs['PseudoTime'] = obs.apply(lambda x : x.PseudoTime1 if np.isnan(x.PseudoTime2) else x.PseudoTime2, axis=1)

    genes = data.index.tolist()
    pool = [[g1, g] for g in genes for g1 in genes]
    adata, TF_names, TF_idx = construct_adata(data, obs, ref_network, time_key='Time', split_k = split_k)
    #adata.obs['cluster'] = cluster['cl'].tolist()
    return adata, ref_network, pool