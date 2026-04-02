#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 22:19:13 2026

@author: chandrima.modak
"""

import os,sys
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse
from pathlib import Path
import warnings
from tqdm import tqdm
import xarray
import pytest

warnings.filterwarnings('ignore')

import preprocess_adata as ppr

def get_sgrna_qc_metrics_macs(crispr_a, min_sgrna_counts=3, q=0.05):
    var_cols = ['sgrna_id','perturbed_gene_name', 'feature_types', 'genome', 'pattern', 'read', 'sequence',
       'n_cells', 'mean_counts', 'total_counts', 'nonz_means']
    # Sanitize excel problems
    crispr_a.var_names = np.where(crispr_a.var_names == '1-Jun', 'JUN-1', crispr_a.var_names)
    crispr_a.var_names = np.where(crispr_a.var_names == '2-Jun', 'JUN-2', crispr_a.var_names)
    # Compute mean of non-zero UMIs
    crispr_a.var['nonz_means'] = ppr._compute_nonzero_means_v1(crispr_a.X)
    sc.pp.calculate_qc_metrics(crispr_a, inplace=True)
    perturb_metadata = crispr_a.var.copy()
    # Annotate perturbed gene
    # perturb_metadata['perturbation_type'] = perturb_metadata.index.str.extract(r'(CRISPRi|CRISPRa|NTC)', expand=False)
    perturb_metadata['perturbed_gene_name'] = perturb_metadata.index.str.replace('-', '_').str.split('_').str[0]
    perturb_metadata['sgrna_id'] = perturb_metadata.index
    perturb_metadata = perturb_metadata.rename(
        {'n_cells_by_counts':'n_cells'}, 
        axis=1) 
    crispr_a.var = perturb_metadata[var_cols].copy()


def process_cellranger_h5_macs(xdata, exp, sample_name, lane):
    '''
    Process single cellranger file
    file_path: h5 or h5ad cellranger filtered data
    exp: experiment name
    sample_name: experiment run
    '''
# =============================================================================
#     file_path = f'{datadir}/CRISPRia_Cellanome_{lane}/per_sample_outs/{sample_name}/count/sample_filtered_feature_bc_matrix.h5'    
#     gex_a, crispr_a = _split_assay(file_path, exp)
#     gex_a.obs['library_id'] = sample_name
#     gex_a.obs['lane_id'] = lane
#     crispr_a.obs['library_id'] = sample_name
#     crispr_a.obs['lane_id'] = lane
#     gex_a.obs_names = gex_a.obs_names + "_" + gex_a.obs['lane_id'] + "_" + gex_a.obs['library_id'] 
#     crispr_a.obs_names = crispr_a.obs_names + "_" + crispr_a.obs['lane_id'] + "_" + crispr_a.obs['library_id']
#     
# =============================================================================
    gex_a, crispr_a = ppr._split_assay(xdata, exp, sample_name, lane)
    # Process sgRNA adata
    get_sgrna_qc_metrics_macs(crispr_a, min_sgrna_counts=3, q=0.05)
    

    # Process scRNA adata
    gex_a.layers['counts'] = gex_a.X.copy()
    # rsc.get.anndata_to_GPU(gex_a)
    gex_a, pre_count, post_count = ppr._basic_qc_gex(gex_a)
    sc.pp.normalize_total(gex_a)
    sc.pp.log1p(gex_a)
    return gex_a, crispr_a, pre_count, post_count