#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate a GRM (nonsparse) and associated phenotypes
Simulate phentotypes off of the ABCD GRM
All of these will be simulated with varying degrees of contributino from the GRM, sites, and noise
Created on Fri Sep 23 10:00:54 2022

@author: christian
"""
import numpy as np
from scipy.stats import random_correlation


#%% Simulate/load a random GRM (nonsparse)

def sim_GRM(n, GRM_save) :
    # Generate a random rorrelation matrix (it needs sum of eigenvalues to equal dimensionality)
    eigs = np.random.rand(n)
    eigs = eigs / sum(eigs) * n
    A = random_correlation.rvs(eigs, tol = 1e-1)
    # Scale correlation into covariance 
    np.save(GRM_save + ".npy", A)


def simulate_phenotypes(GRM, df, sigmas, sim_prefix) :
    num_sites = len(np.unique(df.Zs))
    n = GRM.shape[0]

    
    for i, s in enumerate(sigmas) :
        # Generate random site effects
        Bs = np.random.normal(0,  s[1], num_sites)
        BBs =  {i : Bs[i] for i in range(num_sites)}
        site_effects = [BBs[site] for site in df.Zs]
        
        # compute covariance for simulate(V) for simulation
        V = s[0] * GRM + s[2] * np.identity(n)

        phensuffix = "".join((sigmas[i]*10).astype(int).astype(str))
        print("simulation:  " + phensuffix)
        df.insert(df.shape[1] , sim_prefix + phensuffix, np.random.multivariate_normal(mean = site_effects, cov = V))
    
    return df
    