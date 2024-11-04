import scanpy as sc
from muon import atac as ac


def rna_pre(adata, min_cells=3, target_sum=1e4, n_top_genes=3000):
    adata.raw_data = adata.X
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    return adata


def atac_pre(adata, min_cells=3, target_sum=1e4, n_top_genes=10000):
    adata.raw_data = adata.X
    ac.pp.binarize(adata)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    return adata
