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


def sim_plink_pheno(rng, bed, sigma= [0.5,0 , 0.5], prop_causal = 0.1, npcs=0) :
    nsubjects, nSNPs= bed.shape
    
    # select causal region
    causal_idx = rng.choice(range(nSNPs), size = int(prop_causal * nSNPs), replace=False)
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
    y = y/ np.var(y) * sigma[0]
    
    error = rng.normal(0,1, len(y))
    error = error / np.var(error) *  sigma[2]
    
    y+= error    
    
    return y
    
