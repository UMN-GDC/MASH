#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:48:24 2023

@author: christian
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

#%%


def plot_genetic_clusters(sim, output = None) :
    """
    Plot the clusters from full GRM, shared GRM, non-shared GRM, causal GRM, and non-causal GRM

    Parameters
    ----------
    sim : simulation object
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pca= PCA(n_components=2)
    
    GRM_pcs = pd.DataFrame({"subj_ancestries" : sim.df.subj_ancestries})
    GRM_pcs = GRM_pcs.reset_index()
    
    # full PCs
    GRM_pcs[["PC1_Full", "PC2_Full"]] = pca.fit_transform(sim.genotypes)
    # causal PCs
    GRM_pcs[["PC1_Causal", "PC2_Causal"]] = pca.fit_transform(sim.genotypes[:, sim.causal_idx])
    # noncausal PCs
    GRM_pcs[["PC1_Noncausal", "PC2_Noncausal"]] = pca.fit_transform(sim.genotypes[:, ~sim.causal_idx])
    # shared PCs
    GRM_pcs[["PC1_Shared", "PC2_Shared"]] = pca.fit_transform(sim.genotypes[:, sim.shared_idx])
    # not shared PCs
    GRM_pcs[["PC1_Nonshared", "PC2_Nonshared"]] = pca.fit_transform(sim.genotypes[:, ~sim.shared_idx])
    # shared causal
    GRM_pcs[["PC1_shared_causal", "PC2_shared_causal"]] = pca.fit_transform(sim.genotypes[:, np.intersect1d(sim.shared_idx, sim.causal_idx)])


    test = (pd.wide_to_long(GRM_pcs, stubnames = ["PC1", "PC2"], i = ["index", "subj_ancestries"], 
                            j = "Type", sep = "_", suffix=r'\w+')
            .reset_index())

    p = sns.FacetGrid(test, col="Type", col_wrap = 2)
    p = p.map_dataframe(sns.scatterplot, x = "PC1", y= "PC2", hue = "subj_ancestries")
    
    if output == None :
        p.show()
    else :
        plt.savefig(fname= output + ".png", dpi= 300, format = "png")
