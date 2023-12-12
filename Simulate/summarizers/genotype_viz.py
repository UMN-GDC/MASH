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
# Shared then Non shared

def plotClusters(sim, output = None) :
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
    GRM_pcs[["PC1_Causal", "PC2_Causal"]] = pca.fit_transform(sim.genotypes[:, sim.causals])
    # noncausal PCs
    GRM_pcs[["PC1_Noncausal", "PC2_Noncausal"]] = pca.fit_transform(sim.genotypes[:, ~sim.causals])
    # shared PCs
    shared_idx = np.array(range(int(sim.nSNPs * sim.shared)))
    GRM_pcs[["PC1_Shared", "PC2_Shared"]] = pca.fit_transform(sim.genotypes[:, shared_idx])
    # not shared PCs
    GRM_pcs[["PC1_Nonshared", "PC2_Nonshared"]] = pca.fit_transform(sim.genotypes[:, ~shared_idx])
    # shared causal
    GRM_pcs[["PC1_shared_causal", "PC2_shared_causal"]] = pca.fit_transform(sim.genotypes[:, np.intersect1d(shared_idx, sim.causals)])


    test = (pd.wide_to_long(GRM_pcs, stubnames = ["PC1", "PC2"], i = ["index", "subj_ancestries"], 
                            j = "Type", sep = "_", suffix=r'\w+')
            .reset_index())

    p = sns.FacetGrid(test, col="Type", col_wrap = 2)
    p = p.map_dataframe(sns.scatterplot, x = "PC1", y= "PC2", hue = "subj_ancestries")
    
    if output != None :
        plt.savefig(fname= output + ".png", dpi= 300, format = "png")
    plt.show()



def SNPeffectViz(sim, out =None) :
    """
    Create vizualizations for SNP effects
    Parameters
    ----------
    sim : simulation object
        DESCRIPTION.

    Returns
    -------
    None.
    

    """
    snp_frequencies = np.mean(sim.genotypes[:,sim.causal_idx], axis = 0)/2
    maf = 0.5 - abs(snp_frequencies-0.5)

    # make side by side subplots
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.scatter(maf, abs(sim.SNP_effects))
    ax1.set_xlabel("MAF")
    ax1.set_ylabel("SNP effect size")
    # Seaborn denstiy plot
    sns.kdeplot(ax = ax2, x= abs(sim.SNP_effects))
    ax2.set_xlabel("SNP effect size")
    # increase space between ax1 and ax2 so that the ylabel from ax2 doesn't overlap with ax1
    fig.subplots_adjust(wspace = 0.3)
    if out != None :
        fig.savefig(fname= out, dpi= 300, format = "png")
    fig.show()
    

