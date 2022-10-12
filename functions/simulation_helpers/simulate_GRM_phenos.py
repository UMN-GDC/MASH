#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate a GRM (nonsparse) and associated phenotypes
Simulate phentotypes off of the ABCD GRM
All of these will be simulated with varying degrees of contributino from the GRM, sites, and noise
Created on Fri Sep 23 10:00:54 2022

@author: christian
"""
import os
import numpy as np
from scipy.stats import random_correlation
from scipy.linalg import block_diag


#%% Simulate/load a random GRM (nonsparse)

def sim_GRM(n, GRM_save) :
    if not os.path.exists(GRM_save + ".grm.bin") :
    
        # Generate a random rorrelation matrix (it needs sum of eigenvalues to equal dimensionality)
        eigs = np.random.rand(n)
        eigs = eigs / sum(eigs) * n
        A = random_correlation.rvs(eigs, tol = 1e-1)

        # Create empty upper triangle matrix
        u = np.tril_indices(A.shape[0])

        # select the upper triangle and convert it to float four (standard GCTA GRM format)
        A = A[u].astype('f4')


        # write binary file
        A.tofile(GRM_save + ".grm.bin")

    else: 
        print("No simulations necessary it was already simulated")



def simulate_phenotypes(GRM, df, sigmas, sim_prefix, reps = 1) :
    n = GRM.shape[0]
   
    # Create S similarity matrix 
    sites, sizes= np.unique(df["abcd_site"], return_counts = True)
    #%% Construct the block diagonal
    diags = [np.ones((size,size)) for size in sizes]
    S = np.matrix(block_diag(*diags))

    num_sites = len(sites)
    SG = S * GRM
 
    rng = default_rng()

    for i, s in enumerate(sigmas) :
        # Generate random site effects
        Bs = np.random.normal(0,  s[1], num_sites)
        BBs =  {sites[i] : Bs[i] for i in range(num_sites)}
        site_effects = [BBs[site] for site in df.abcd_site]

        # compute covariance for simulate(V) for simulation
        V = s[0] * GRM + s[2] * SG +  s[3] * np.identity(n)

        paramset = "".join((sigmas[i]*10).astype(int).astype(str))

	# repeat it for the number of desired replicates
        rands = pd.DataFrame(rng.multivariate_normal(mean = site_effects, cov = V, reps).T)
        rands.columns = [paramset + "_" + str(i+1) for i in range(reps)] 
    
    df = pd.concat([df, rands], axis = 1)
    
    return df
    
