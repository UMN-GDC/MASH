#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:48:24 2023

@author: christian
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def get_cluster_pcs(sim):
    """
    Get PCA components for different genotype subsets

    Parameters
    ----------
    sim : simulation object

    Returns
    -------
    pd.DataFrame with PC columns for different genotype types
    """
    nsubjects, nSNPs = sim.genotypes.shape
    pca = PCA(n_components=2)
    
    GRM_pcs = pd.DataFrame({"subj_ancestries": sim.df.subj_ancestries})
    GRM_pcs = GRM_pcs.reset_index()
    
    GRM_pcs[["PC1_Full", "PC2_Full"]] = pca.fit_transform(sim.genotypes)
    
    return GRM_pcs

