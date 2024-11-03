from collections import defaultdict

import scanpy as sc
from sklearn.neighbors import kneighbors_graph
from sklearn import decomposition
import numpy as np
from scipy.spatial import distance
import muon as mu

from scmfg.utils import rna_pre, atac_pre


class SCMFG:
    def __init__(
            self,
            mdata,
            n_components: int = 25,
            latent_dim: int = 20,
            mofa_factor: int = 10,
            is_preprocessed: bool = False,
    ):
        """Initialize a SCMFG object for analysis

        Args:
            mdata (md.MuData):
                The input MuData object.
            n_components (int, optional)
                The number of feature groups per omics. Default to 15.
            latent_dim (int, optional):
                The latent dimension of the model. Defaults to 20.
            mofa_factor (int, optional)
                The number of MOFA factors. Default to 10.
            is_preprocessed (bool, optional):
                If True, skip data preprocessing operations.
        """
        self.mdata = mdata
        self.mods = list(mdata.mod.keys())
        self.n_components = n_components
        self.latent_dim = latent_dim
        self.mofa_factor = mofa_factor
        self.is_preprocessed = is_preprocessed

    def run(self):

        # Preprocess the input Mudata data
        mdata_pre = self.preprocess()

        # Each omics is grouped using the LDA model to get the group category to which the feature belongs
        label = {}
        for mod in self.mods:
            label[mod] = self.get_label(mdata_pre[mod])

        # Calculate the K-nearest neighbors of each cell based on each feature group
        knn = {}
        for mod in self.mods:
            knn[mod] = self.get_knn(mdata_pre[mod], label[mod])

        # Calculate the correlation between different omics feature groups
        knn_corr = self.get_knn_corr(knn)

        # Integrate the most similar feature groups from different omics
        self.mdata.obsm["X_scMFG"] = self.integrate(mdata_pre, label, knn_corr)

    def preprocess(self):
        mdata_pre = {}
        for mod in self.mods:
            adata = self.mdata[mod].copy()
            adata.var_names_make_unique()
            if self.is_preprocessed:
                mdata_pre[mod] = adata
                continue

            if mod == self.mods[0]:
                # Preprocessed scRNA data
                adata = rna_pre(adata)
            elif mod == self.mods[1]:
                # Preprocessed scATAC data
                adata = atac_pre(adata)
            # Select highly variable features
            mdata_pre[mod] = adata[:, adata.var["highly_variable"]]
        return mdata_pre

    def get_label(self, features):

        # Each omics is grouped using the LDA model
        lda_model = decomposition.LatentDirichletAllocation(n_components=self.n_components)
        topic = lda_model.fit_transform(features.X.T)
        group = np.zeros((features.shape[1], 1), dtype=int)

        # Find the group to which the feature belongs
        for i in range(features.X.shape[1]):
            group[i] = np.argmax(topic[i, :])
        label = dict([(k, []) for k in range(self.n_components)])
        for i in range(self.n_components):
            for j in range(group.shape[0]):
                if group[j] == i:
                    label[i].append(j)
        return label

    def get_knn(self, adata, label, k_num=15):
        knn = {k: [] for k in ["group_" + str(i) for i in range(self.n_components)]}
        for i in range(self.n_components):
            if label[i]:
                # Compute KNN of cells for a specific feature group
                knn_dist = kneighbors_graph(adata.X[:, label[i]], mode="connectivity", n_neighbors=k_num)
                knn_dist_coo = knn_dist.tocoo()
                neighbors_dict = defaultdict(list)
                for j, k in zip(knn_dist_coo.row, knn_dist_coo.col):
                    neighbors_dict[j].append(k)

                for j in neighbors_dict:
                    knn["group_" + str(i)].extend(neighbors_dict[j][:k_num])
        return knn

    def get_knn_corr(self, knn):
        knn_corr = {k: 0 for k in ["group_" + str(i) for i in range(self.n_components)]}
        tmp = []
        # The similarity between feature groups was calculated using jaccard
        for i in range(self.n_components):
            for j in range(self.n_components):
                tmp.append(1 - distance.jaccard(knn[self.mods[0]]["group_" + str(i)],
                                                knn[self.mods[1]]["group_" + str(j)]))
            max_k = 0
            flag = 0
            # Find the most similar feature groups
            for k in range(self.n_components):
                if tmp[k] > max_k and k not in knn_corr.values():
                    max_k = tmp[k]
                    flag = k
            knn_corr["group_" + str(i)] = flag
            tmp = []
        return knn_corr

    def integrate(self, mdata, label, knn_corr):
        var_names = []
        rna_name = self.mods[0]
        atac_name = self.mods[1]
        for i in range(self.n_components):
            if label[rna_name][i] and label[atac_name][i]:
                var_name = mu.MuData({
                    rna_name: mdata[rna_name][:, label[rna_name][i]],
                    atac_name: mdata[atac_name][:, label[atac_name][knn_corr["group_" + str(i)]]]})
                if var_name[i]:
                    # Use mofa to integrate feature groups between different omics
                    mu.tl.mofa(var_name, n_factors=self.mofa_factor)
                var_names.append(var_name)
        # Save the mofa integration results
        integrate_dim = []
        for i in range(self.n_components):
            if label[rna_name][i] and label[atac_name][i]:
                for j in range(self.mofa_factor):
                    integrate_dim.append(var_names[i].obsm["X_mofa"][:, j])
        integrate_array = np.asarray(integrate_dim).T
        # Dimensionality reduction using pca
        return sc.tl.pca(integrate_array, n_comps=self.latent_dim)
