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
import pandas as pd
from scipy.stats import random_correlation
from scipy.linalg import block_diag
from functions.load_data import ReadGRMBin
from functions.AdjHE_estimator import load_n_estimate



#%% Simulate/load a random GRM (nonsparse)

def sim_GRM(n) :
    #if not os.path.exists(GRM_save + ".grm.bin") :
    
    # Generate a random rorrelation matrix (it needs sum of eigenvalues to equal dimensionality)
    eigs = np.random.rand(n)
    eigs = eigs / sum(eigs) * n
    A = random_correlation.rvs(eigs, tol = 1e-3)
    # Standardize it 
    #A = n/1.5 * A/ A.sum()

    return A

        
def sim_pheno(n, sigma) :
    # make sure no zero variance components
    for i, v in enumerate(sigma) :
        if v == 0:
            sigma[i] = 1e-6

    df = pd.DataFrame({"fid" : np.arange(n),
                       "iid" : np.arange(n)})
    
    rng = np.random.default_rng()
    df["abcd_site"] = rng.choice(np.arange(25), size= (n,))
    
    GRM = sim_GRM(n)
 
    # Create S similarity matrix 
    sites, sizes= np.unique(df["abcd_site"], return_counts = True)
    # Construct the block diagonal
    diags = [np.ones((size,size)) for size in sizes]
    S = block_diag(*diags)
    # Standardize S
    #S = n*S/ S.sum()

    # Create interaction matrix
    SG = S * GRM
    # Standardize it
    #SG = SG/ SG.sum()

    I = np.eye(n)    

    df["G"] = rng.multivariate_normal(mean = np.repeat(0,n), cov = GRM)    
    df["S"] = rng.multivariate_normal(mean = np.repeat(0,n), cov = S)    
    df["SG"] = rng.multivariate_normal(mean = np.repeat(0,n), cov = SG)    
    df["I"] = rng.multivariate_normal(mean = np.repeat(0,n), cov =  I)


    # Calculate desired variance ratios
    StoG = sigma[1]/sigma[0]
    SGtoG = sigma[2]/sigma[0]
    EtoG = sigma[3]/sigma[0]


    # Calculate empricial variance ratios
    gen_var = np.var(df["G"])
    site_var = np.var(df["S"])
    sg_var = np.var(df["SG"])
    e_var = np.var(df["I"])


    # Find the empirical ratio of variances
    StoG_sim = site_var / gen_var
    site_variance_scaling = StoG_sim / StoG
    
    SGtoG_sim = sg_var / gen_var
    sg_variance_scaling = SGtoG_sim / SGtoG
    
    EtoG_sim = e_var / gen_var
    error_variance_scaling = EtoG_sim / EtoG

    # Scale site effects so total contributed variance is as presecribed
    df["S"] = df["S"] / np.sqrt(site_variance_scaling)
    df["SG"] = df["SG"] / np.sqrt(sg_variance_scaling)
    df["I"] = df["I"] / np.sqrt(error_variance_scaling)

    # Scale by their contributions    
    df["Y"] = df.G + df.S + df.SG + df.I
    
    return GRM, df


#%%

GRM, df = sim_pheno(500, sigma = [0.5, 0.125, 0.125, 0.25])
load_n_estimate(df=df, covars= [], nnpc=0, mp = "Y", GRM = GRM, std = False, fast= True, RV = None)


#%%
#%% Simulations
import itertools
#import plotly.express as px
#from plotly.offline import plot

sg = [0.25, 0.5, 0.75]
ss = [0, 0.25]
ssg = [0, 0.25]
se = [0, 0.25, 0.5, 0.75]

sigmas = []
for ses in itertools.product(sg,ss, ssg, se) :
    if (sum(ses) == 1) and (ses[0] > 0.1) :
        sigmas.append(ses)
      
# Add a completely null case
sigmas.insert(0, (0, 0,0, 1))

# remove nonsesne setups
sigmas = np.array(sigmas)[[0, 1,3,4,5, 7, 9]]

#%%
results = {"sg" : [], 
           "ss" : [],
           "ssg" : [],
           "se" : [],
           "Basic_est" : [],
           "Site_RE" : [],
           "Site_FE" : []}
for sigma in sigmas :
    for rep in range(20) :
        results["sg"].append(sigma[0])
        results["ss"].append(sigma[1])
        results["ssg"].append(sigma[2])
        results["se"].append(sigma[3])
   
        
        GRM, df = sim_pheno(n= 100, sigma = list(sigma))
        # Fit RE AdjHE
        results["Site_RE"].append(load_n_estimate(df=df, covars= [], nnpc=0, mp = "Y",
                                          GRM = GRM, std = False, fast= True, RV = "abcd_site")["h2"][0])
    
        # Fit basic AdjHE w/o site
        results["Basic_est"].append(load_n_estimate(df=df, covars= [], nnpc=0, mp = "Y",
                                          GRM = GRM, std = False, fast= True, RV = None)["h2"][0])
        
        # Fit AdjHE with Site fixed effect
        results["Site_FE"].append(load_n_estimate(df=df, covars= ["abcd_site"], nnpc=0, mp = "Y",
                                          GRM = GRM, std = False, fast= True, RV = None)["h2"][0])

        
results = pd.DataFrame(results)
results2 = pd.melt(results, id_vars= ["sg", "ss", "se"], value_vars=['Basic_est', 'Site_RE', 'Site_FE'])
results2["herit"] = results2.sg/ (results2.sg + results2.ss + results2.se)
#%%
results.to_csv("sim_1000_10s.csv",index = False)



