from scipy.stats import spearmanr
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, r2_score, roc_auc_score
from sklearn.neighbors import KNeighborsRegressor
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from scipy.stats import norm
import torchdiffeq


def pred_density_fate_bias(x0_adata, adata, time_key, embedding_key, cell_type_key, model, target_cell_types):
    pred_x = np.zeros_like(x0_adata.obsm[embedding_key])
    t_idxs = [x0_adata.obs[time_key] == t for t in np.sort(np.unique(x0_adata.obs[time_key]))]
    for j, t in enumerate(np.sort(np.unique(x0_adata.obs[time_key]))):
        pred_x[t_idxs[j]] = torchdiffeq.odeint(model, torch.Tensor(x0_adata.obsm[embedding_key][t_idxs[j]]), t=torch.Tensor([t,2]), method='rk4')[-1].detach().numpy()
    fate_bias = get_density_fate_bias(pred_x, adata, embedding_key, cell_type_key, target_cell_types)
    return fate_bias

def get_density_fate_bias(pred_x, adata, embedding_key, cell_type_key, target_cell_types):
    mono_mean = np.mean(adata[adata.obs[cell_type_key]==target_cell_types[1]].obsm[embedding_key], axis=0)
    mono_std = np.std(adata[adata.obs[cell_type_key]==target_cell_types[1]].obsm[embedding_key], axis=0)

    neu_mean = np.mean(adata[adata.obs[cell_type_key]==target_cell_types[0]].obsm[embedding_key], axis=0)
    neu_std = np.std(adata[adata.obs[cell_type_key]==target_cell_types[0]].obsm[embedding_key], axis=0)

    mono_fate = np.sum(np.log(norm(mono_mean, mono_std).pdf(pred_x)), axis=1)
    neu_fate = np.sum(np.log(norm(neu_mean,  neu_std).pdf(pred_x)), axis=1)
    fate_bias = np.exp(neu_fate) / (np.exp(neu_fate) + np.exp(mono_fate) )
    return fate_bias


def fate_map(adata, bias_list = ['Reprogrammed', 'Failed'] ):
    adata.obs['got_fate_bias'] = adata.obsm['descendant'][bias_list[0]] / (adata.obsm['descendant'][bias_list[0]] + adata.obsm['descendant'][bias_list[1]])
    adata.obs.loc[np.isnan(adata.obs.got_fate_bias), 'got_fate_bias'] = 0.5
    
    return adata.obs['got_fate_bias']


def get_fate_bias(pred_cell_fate, bias_list = ['Reprogrammed', 'Failed']):
    fate_bias = (pred_cell_fate[bias_list[0]].to_numpy() / (pred_cell_fate[bias_list[0]].to_numpy() + pred_cell_fate[bias_list[1]]).to_numpy())
    fate_bias[np.isnan(fate_bias)] = .5
    return fate_bias

def smoothe_score(X, y, n_neighbors=5):
    idx = ~np.isnan(y)
    knn = KNeighborsRegressor(n_neighbors=n_neighbors)
    knn.fit(X[idx], y[idx])
    y_smoothed = knn.predict(X)
    return y_smoothed, knn


def correlation_2x2(pred_fate_bias, x0_adata, fate_bias_key, embedding_key):
    pred_fate_bias[np.isnan(pred_fate_bias)] = .5
    x0_adata.obs['pred_fate_bias'] = pred_fate_bias
    x0_adata.obs['smoothed_fate_bias'], _ = smoothe_score(x0_adata.obsm[embedding_key], x0_adata.obs[fate_bias_key])
    x0_adata.obs['smoothed_pred_fate_bias'], _ = smoothe_score(x0_adata.obsm[embedding_key], x0_adata.obs['pred_fate_bias'])
    spearman = []
    
    spearman.append(spearmanr(x0_adata.obs['pred_fate_bias'], x0_adata.obs[fate_bias_key])[0])
    spearman.append(spearmanr(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs[fate_bias_key])[0])
    spearman.append(spearmanr(x0_adata.obs['pred_fate_bias'], x0_adata.obs['smoothed_fate_bias'])[0])
    spearman.append(spearmanr(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs['smoothed_fate_bias'])[0])
    
    pearson = []
    pearson.append(pearsonr(x0_adata.obs['pred_fate_bias'], x0_adata.obs[fate_bias_key])[0])
    pearson.append(pearsonr(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs[fate_bias_key])[0])
    pearson.append(pearsonr(x0_adata.obs['pred_fate_bias'], x0_adata.obs['smoothed_fate_bias'])[0])
    pearson.append(pearsonr(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs['smoothed_fate_bias'])[0])
    
    mse = []
    mse.append(mean_squared_error(x0_adata.obs['pred_fate_bias'], x0_adata.obs[fate_bias_key]))
    mse.append(mean_squared_error(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs[fate_bias_key]))
    mse.append(mean_squared_error(x0_adata.obs['pred_fate_bias'], x0_adata.obs['smoothed_fate_bias']))
    mse.append(mean_squared_error(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs['smoothed_fate_bias']))
    
    
    r2 = []
    r2.append(r2_score(x0_adata.obs['pred_fate_bias'] , x0_adata.obs[fate_bias_key]))
    r2.append(r2_score(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs[fate_bias_key]))
    r2.append(r2_score(x0_adata.obs['pred_fate_bias'], x0_adata.obs['smoothed_fate_bias']))
    r2.append(r2_score(x0_adata.obs['smoothed_pred_fate_bias'], x0_adata.obs['smoothed_fate_bias']))
    
    
    accs = []
    accs.append(roc_auc_score(x0_adata.obs[fate_bias_key] > 0.5, x0_adata.obs['pred_fate_bias'] ))
    accs.append(roc_auc_score(x0_adata.obs[fate_bias_key] > 0.5, x0_adata.obs['smoothed_pred_fate_bias']))
    accs.append(roc_auc_score(x0_adata.obs['smoothed_fate_bias'] > 0.5, x0_adata.obs['pred_fate_bias']))
    accs.append(roc_auc_score(x0_adata.obs['smoothed_fate_bias'] > 0.5, x0_adata.obs['smoothed_pred_fate_bias']))

    return pd.DataFrame([spearman, pearson, mse, r2, accs], index=['scc', 'pcc', 'mse', 'r2', 'auc'], columns=['pred vs gt', 'sm pred vs gt', 'pred vs sm gt', 'sm pred vs sm gt'])

def positive_control(x0_adata, embedding_key, groudtruth_key, n_splits=5):
    X = x0_adata.obsm[embedding_key][:, ]
    linear_model = LinearRegression()

    y = x0_adata.obs[groudtruth_key].values
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=42)  # 5折交叉验证
    y_pred_cv = cross_val_predict(linear_model, X, y, cv=cv)


    mse = mean_squared_error(y, y_pred_cv)
    r2 = r2_score(y, y_pred_cv)
    auc = roc_auc_score(y > 0.5, y_pred_cv)
    pcc =  pearsonr(y, y_pred_cv)[0]
    scc = spearmanr(y, y_pred_cv)[0]
    print("SCC:", scc)
    print("PCC:", pcc)
    print("Mean Squared Error:", mse)
    print("R-squared:", r2)
    print("ROC-AUC:", auc)
    return scc, pcc, mse, r2, auc