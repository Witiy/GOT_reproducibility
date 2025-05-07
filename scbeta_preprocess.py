import logging
import numpy as np
import scipy as sp
from typing import *


import pandas as pd

import igraph as ig
from sklearn import preprocessing


class CountsNormalizer:
    """
    Normalize and optionally standardize a dataset, dealing properly with edge cases such as division by zero.
    """
    def __init__(self,
            mu: np.ndarray = None,
            sd: np.ndarray = None,
            totals = None,
            mean_center: bool = True,
            normalize_variance: bool = True,
            level: int = 10000) -> None:

        self.sd = sd  # type: np.ndarray
        self.mu = mu  # type: np.ndarray
        self.totals = totals  # type: np.ndarray
        self.level = level
        self.mean_center = mean_center
        self.normalize_variance = normalize_variance

    def norm_factor_from_totals(self, totals):
        norm_factor = (self.level/totals)
        norm_factor[~np.isfinite(norm_factor)] = 0
        norm_factor = sp.sparse.diags(norm_factor).tocsr()
        return norm_factor

    def fit(self, vals: sp.sparse.csr_matrix, cells: np.ndarray = None) -> None:
        totals = vals.sum(0).A.reshape(-1)
        norm_factor = self.norm_factor_from_totals(totals)

        vals = vals.dot(norm_factor).T
        if cells is not None:
            vals = vals[cells, :]

        #Standard scaler is particularly fast and works with sparse arrays (it computes the mean even if it is not used in rescaling).
        scaler = preprocessing.StandardScaler(with_mean=False)
        scaler.fit(vals)

        self.mu = scaler.mean_
        self.sd = np.sqrt(scaler.var_)
        self.totals = totals

    def transform(self, vals: sp.sparse.csr_matrix, cells: np.ndarray = None, genes: np.ndarray = None) -> np.ndarray:
        """
        Normalize a matrix using the previously calculated aggregate statistics

        Args:
            vals (sp.sparse.csr_matrix):     Matrix of shape (n_genes, n_cells)
            cells (ndarray):    Optional indices of the cells that are represented in vals
            cells (ndarray):    Optional indices of the genes that are represented in vals

        Returns:
            vals_adjusted (ndarray):    The normalized values
        """

        if genes is not None:
            vals = vals[genes, :]
        else:
            genes = np.arange(vals.shape[0])

        norm_factor = norm_factor = self.norm_factor_from_totals(self.totals)
        if cells is not None:
            vals = vals[:, cells]
            norm_factor = norm_factor[cells, :][:, cells]

        vals = vals.dot(norm_factor).A

        if self.mean_center:
            vals = vals - self.mu[genes][:, None]
        if self.normalize_variance:
            vals = div0(vals.T, self.sd[genes]).T
            
        return vals
    
def div0(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


from sklearn.svm import SVR

class FeatureSelection:
    def __init__(self) -> None:
        self.genes = None  # type: np.ndarray
        self.mu = None  # type: np.ndarray
        self.sd = None  # type: np.ndarray
        self.totals = None  # type: np.ndarray

    def fit(self,  vals: sp.sparse.csr_matrix,
            cells: np.ndarray = None,
            mu: np.ndarray = None,
            sd: np.ndarray = None) -> np.ndarray:
        """
        Fits a noise model (CV vs mean)

        Args:
            ds (LoomConnection):    Dataset
            n_genes (int):  number of genes to include
            cells (ndarray): cells to include when computing mean and CV (or None)
            mu, std:        Precomputed mean and standard deviations (optional)

        Returns:
            ndarray of selected genes (list of ints)
        """
        if cells is not None:
            vals = vals[:, cells]
        self.fit_cells = cells

        if mu is None or sd is None:
            scaler = preprocessing.StandardScaler(with_mean=False)
            scaler.fit(vals.T)
            mu = scaler.mean_
            sd = np.sqrt(scaler.var_)

        # if "_Valid" in ds.ra:
        #     valid = ds.ra._Valid == 1
        # else:
        #     valid = np.ones(ds.shape[0], dtype='bool')
        # valid = ((valid) & (mu > self.min_expr)).astype(int)
        ok = (mu > 0) & (sd > 0)

        self.mu = mu
        self.sd = sd
        self.ok = ok
        # self.valid = valid

        cv = sd[ok] / mu[ok]
        log2_m = np.log2(mu[ok])
        log2_cv = np.log2(cv)

        svr_gamma = 1000. / len(mu[ok])
        clf = SVR(gamma=svr_gamma)
        clf.fit(log2_m[:, np.newaxis], log2_cv)
        fitted_fun = clf.predict
        # Score is the relative position with respect of the fitted curve
        fitted_val = fitted_fun(log2_m[:, np.newaxis])
        score = log2_cv - fitted_val
        # score = score * valid[ok]

        self.fitted_val = fitted_val
        self.score = score
        self.log2_cv = log2_cv
        self.log2_m = log2_m

    def select_genes(self, n_genes=1000, min_expr=0.001, valid_genes=None):
        _valid = (self.mu > min_expr)
        if valid_genes is not None:
            _valid = _valid & valid_genes
        _valid_score = self.score * _valid[self.ok].astype(float)

        picked_genes = np.where(self.ok)[0][np.argsort(_valid_score)][-n_genes: ]
        return picked_genes


from sklearn.decomposition import PCA

class PCAProjection:
    """
    Project a dataset into a reduced feature space using PCA. The projection can be fit
    to one dataset then used to project another. To work properly, both datasets must be normalized in the same 
    way prior to projection.
    """
    def __init__(self, genes: np.ndarray, max_n_components: int = 50) -> None:
        """
        Args:
            genes:              The genes to use for the projection
            max_n_components:   The maximum number of projected components
            nng                 Non-neuronal genes, to be zeroed in neurons (where TaxonomyRank1 == "Neurons")
        """
        self.genes = genes
        self.n_components = max_n_components

        self.cells = None  # type: np.ndarray
        self.pca = None  # type: IncrementalPCA
        self.sigs = None  # type: np.ndarray
        # self.scan_batch_size = scan_batch_size

    def fit(self, vals: sp.sparse.csr_matrix, normalizer, cells: np.ndarray = None) -> None:
        n_cells = vals.shape[1] if cells is None else cells.shape[0]
        n_genes = self.genes.shape[0]

        self.pca = PCA(n_components=self.n_components)
        norm_vals = normalizer.transform(vals, genes=self.genes, cells=cells)
        self.pca.fit(norm_vals.T)

    def transform(self, vals: sp.sparse.csr_matrix, normalizer, cells: np.ndarray = None) -> np.ndarray:
        
        n_cells = vals.shape[1] if cells is None else cells.shape[0]
        n_genes = self.genes.shape[0]

        norm_vals = normalizer.transform(vals, genes=self.genes, cells=cells)
        transformed = self.pca.transform(norm_vals.T)
        return transformed

    def fit_transform(self, vals: sp.sparse.csr_matrix, normalizer, cells: np.ndarray = None) -> np.ndarray:
        self.fit(vals, normalizer, cells)
        return self.transform(vals, normalizer, cells)