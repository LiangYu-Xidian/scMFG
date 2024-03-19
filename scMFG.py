# import hotspot
import sys

sys.path.append('C:/Users/21593/PycharmProjects/scmer/sourcecode')
from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import spectral_clustering
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from scipy.sparse import issparse
from mudata import MuData
from sklearn.decomposition import TruncatedSVD
from sklearn import decomposition
import scanpy as sc
import numpy as np
# from scMVP.dataset import GeneExpressionDataset, CellMeasurement
from scipy.spatial import distance
# import cosg as cosg
import pandas as pd
import anndata as ad
import networkx as nx
import random
from scipy.sparse import csr_matrix
import os
import math
import scipy.stats as st
import copy
import muon as mu
from muon import atac as ac
from sklearn.mixture import GaussianMixture as GMM
# import community as community_louvain
from matplotlib import RcParams
import matplotlib.pyplot as plt
from scipy import stats
from sklearn import metrics

np.random.seed(0)


# from hyperopt import hp


def readData():
    print(os.getcwd())
    rna_file = "./data/data_h5ad/share_skin/RNA_unpred4.h5ad"
    # rna_file = "./data/data_h5ad/share_skin_hvg/RNA_unpred4.h5ad"
    rna = sc.read_h5ad(rna_file)
    dataset = "share_skin"
    atac_file = "./data/data_h5ad/share_skin/ATAC_unpred.h5ad"
    # atac_file = "./data/data_h5ad/share_skin_hvg/ATAC_unpred.h5ad"
    atac = sc.read_h5ad(atac_file)
    return rna, atac, dataset


def rna_pre(adata):
    adata.raw_data = adata.X
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    # scran = adata.X / adata.obs["size_factors"].values[:, None]
    #
    # adata.X = csr_matrix(sc.pp.log1p(scran))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=5000)
    return adata


def atac_pre(adata):
    adata.raw_data = adata.X
    ac.pp.binarize(adata)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=10000)
    return adata


if __name__ == "__main__":

    rna, atac, dataname = readData()
    print("################")
    print(dataname)
    celltype = {"sciCAR_cellline": "labels", "snare_p0": "cell_type", "snare_cellline": "cell_line",
                "share_skin": "celltype", "10x_pbmc": "label", "kidney": "cell_name", "neuips": "cell_type",
                "10x_lymph_node": "cell_type", "snare_AdBrainCortex": "cell_type"}

    cell_label = celltype[dataname]
    rna.var_names_make_unique()
    atac.var_names_make_unique()
    rna = rna_pre(rna)
    atac = atac_pre(atac)
    feature_rna = rna[:, rna.var["highly_variable"]]
    feature_atac = atac[:, atac.var["highly_variable"]]
    print("feature_rna:", feature_rna.X.shape)
    print("feature_atac:", feature_atac.X.shape)
    comp = 25
    print("comp:", comp)

    lda_model_rna = decomposition.LatentDirichletAllocation(n_components=comp)
    rna_topic = lda_model_rna.fit_transform(feature_rna.X.T)


    lda_model_atac = decomposition.LatentDirichletAllocation(n_components=comp)
    atac_topic = lda_model_atac.fit_transform(feature_atac.X.T)


    group_rna = np.zeros((feature_rna.shape[1], 1), dtype=int)
    for i in range(feature_rna.X.shape[1]):
        group_rna[i] = np.argmax(rna_topic[i, :])

    group_atac = np.zeros((feature_atac.shape[1], 1), dtype=int)
    for i in range(feature_atac.X.shape[1]):
        group_atac[i] = np.argmax(atac_topic[i, :])
        key = []
    for i in range(comp):
        key.append(i)


    rna_label = dict([(k, []) for k in key])
    for i in range(comp):
        for j in range(group_rna.shape[0]):
            if group_rna[j] == i:
                rna_label[i].append(j)
        key = []
    for i in range(comp):
        key.append(i)
    atac_label = dict([(k, []) for k in key])
    for i in range(comp):
        for j in range(group_atac.shape[0]):
            if group_atac[j] == i:
                atac_label[i].append(j)

    k_num = 15
    print("k_num:", k_num)
    # from sklearn.neighbors import kneighbors_graph
    # key = []
    # for i in range(rna.X.shape[0]):
    #     key.append("cell_" + str(i))
    # cell_rna = dict([(k, []) for k in key])
    key = []
    for i in range(comp):
        key.append("group_" + str(i))
    knn_rna = dict([(k, []) for k in key])
    for i in range(comp):
        if rna_label[i] == []:
            pass
        else:
            rna_knn_dist = kneighbors_graph(rna.X[:, rna_label[i]], mode="connectivity", n_neighbors=k_num)
            rna_knn_dist_array = rna_knn_dist.toarray()
            for j in range(rna.X.shape[0]):
                for k in range(rna.X.shape[0]):
                    if rna_knn_dist_array[j][k] == 1:
                        knn_rna["group_" + str(i)].append(k)
                # print(count)
            # print(len(knn_rna["group_" + str(i)]))

    key = []
    for i in range(comp):
        key.append("group_" + str(i))
    knn_atac = dict([(k, []) for k in key])
    for i in range(comp):
        if atac_label[i] == []:
            pass
        else:
            atac_knn_dist = kneighbors_graph(atac.X[:, atac_label[i]], mode="connectivity", n_neighbors=k_num)
            atac_knn_dist_array = atac_knn_dist.toarray()
            for j in range(atac.X.shape[0]):
                for k in range(atac.X.shape[0]):
                    if atac_knn_dist_array[j][k] == 1:
                        knn_atac["group_" + str(i)].append(k)


    key = []
    for i in range(comp):
        key.append("group_" + str(i))
    knn_corre = dict([(k, 0) for k in key])
    tmp = []
    for i in range(comp):
        for j in range(comp):
            tmp.append(1 - distance.jaccard(knn_rna["group_" + str(i)], knn_atac["group_" + str(j)]))
        # print(tmp)
        max = 0
        flag = 0
        for k in range(comp):
            if tmp[k] > max and k not in knn_corre.values():
                max = tmp[k]
                flag = k
            # print(distance.jaccard(knn_rna["group_" + str(i)],knn_atac["group_" + str(j)]))
        knn_corre["group_" + str(i)] = flag
        tmp = []


    var_name = [None for i in range(comp)]
    for i in range(comp):
        # print(i)
        if rna_label[i] == [] or atac_label[i] == []:
            pass
        else:
            var_name[i] = mu.MuData({"rna": feature_rna[:, rna_label[i]],
                                     "atac": feature_atac[:, atac_label[knn_corre["group_" + str(i)]]]})
            mu.tl.mofa(var_name[i], n_factors=factor)
    dim_inte = []
    for i in range(comp):
        if rna_label[i] == [] or atac_label[i] == []:
            pass
        else:
            for j in range(factor):
                dim_inte.append(var_name[i].obsm["X_mofa"][:, j])
    array_inte = np.asarray(dim_inte).T
    array_inte.shape
    pca_dim = 20
    # print(array_inte)
    pca_inte = sc.tl.pca(array_inte, n_comps=pca_dim)
    rna.obsm['scMFG'] = pca_inte


