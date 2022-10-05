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


#%% Simulate/load a random GRM (nonsparse)

def sim_GRM(n, GRM_save) :
    if not os.path.exists(GRM_save) :
    
        # Generate a random rorrelation matrix (it needs sum of eigenvalues to equal dimensionality)
        eigs = np.random.rand(n)
        eigs = eigs / sum(eigs) * n
        A = random_correlation.rvs(eigs, tol = 1e-1)
        # Scale correlation into covariance 
        np.save(GRM_save + ".npy", A)
    else: 
        print("No simulations necessary it was already simulated")

def simulate_phenotypes(GRM, df, sigmas, sim_prefix, reps = 1) :
    sites = np.unique(df.abcd_site)
    n = GRM.shape[0]
    num_sites = len(sites)
    
    for i, s in enumerate(sigmas) :
        # Generate random site effects
        Bs = np.random.normal(0,  s[1], num_sites)
        BBs =  {sites[i] : Bs[i] for i in range(num_sites)}
        site_effects = [BBs[site] for site in df.abcd_site]
        
        # compute covariance for simulate(V) for simulation
        V = s[0] * GRM + s[2] * np.identity(n)

        paramset = "".join((sigmas[i]*10).astype(int).astype(str))

	# repeat it for the number of desired replicates
        for i in range(reps): 
            colname = paramset + "_" + str(i+1)
            print("simulation:  " + colname)
            df.insert(df.shape[1] , sim_prefix + colname, np.random.multivariate_normal(mean = site_effects, cov = V))
    
    return df
    
