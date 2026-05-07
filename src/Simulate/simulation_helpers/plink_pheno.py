#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 13:17:52 2023

@author: christian
"""
import numpy as np
from statsmodels.regression.linear_model import OLS
from sklearn.decomposition import PCA
import dask as di


rng = np.random.default_rng()


def sim_plink_pheno(rng, bed, h2= 0.5, prop_causal = 0.1, npcs=0, causals = None) :
    """
    Simulate a phenotype from a plink files.

    Parameters
    ----------
    rng : random number generatorng
        numpy random number generator.
    bed: array
        ( x ) array of unstandardized genotypes read in from bed file.
    h2 : float, optional
        heritability of the phenotype. The default is 0.5.
    causals : list, optional
        list of causal snp positions indexed to the genotype array (i.e. it's not absolute position in the genomic sense
        but rather the position along the bim file.
    prop_causal : float, optional
        if causals are unspecified, then specify the number of SNPs to randomly select as causal. The default is 0.1.
    npcs : int, optional
        number of PCs to subtract from the phenotype. The default is 0.
    Returns
    -------
    Gene_contrib : array 
        1-D array containing the contribution of genetics to each subjects phenotype.
    causals : array
        1-D array of the genotype index that was selected to be causal
    snp_effects: array
        1-D array of SNP effects
    """
    nsubjects, nSNPs= bed.shape
    
    # select causal region if not specified
    if causals is None : 
        causal_idx = rng.choice(range(nSNPs), size = int(prop_causal * nSNPs), replace=False)
    else :
        causal_idx = causals

    noncausal_idx = []
    for i in range(nSNPs): 
        if i not in causal_idx : 
            noncausal_idx.append(i)
            
    # Simulate betas
    # beta_non= 0 # place holder in case we wantto simulate it from different distribution
    beta_causal = rng.normal(3,1, size = len(causal_idx))
    
    betas = np.zeros((nSNPs, 1))
    for i, c_index in enumerate(causal_idx):
        betas[c_index,0] = beta_causal[i]
    
    # create phenotype
    y = np.dot(bed, betas)
    # Calculate PC's
    pcs = PCA(n_components = 5).fit_transform(bed)
    
    # subtract contributions of the PC's
    y = np.array(y)    

    y = OLS(endog = y, exog= pcs).fit().resid
    y = y/ np.var(y) * h2 
    
    error = rng.normal(0,1, len(y))
    error = error / np.var(error) * (1-h2) 
    
    y+= error    
    
    return y
    
