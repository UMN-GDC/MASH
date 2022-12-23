#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 16:23:04 2022

@author: christian
"""
from functions.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pandas as pd
from functions.Estimation.all_estimators import Basu_estimation
from run_sims import sim_experiment, sim_n_est
import matplotlib.pyplot as plt
import seaborn as sns

#%%



def sim_n_est(nsubjects = 1000, sigma = [0.5,0.25, 0.25], site_comp = "IID", nsites = 30,
              theta_alleles =0.5, nclusts =1, dominance=3, prop_causal=0.25, site_dep=False, nnpc = 1,
              nSNPs=20, phens = 2, site_het = False) :
    sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
    # Run through full simulation and estimation
    sim.full_sim(nsites= nsites, sigma= sigma, phens = phens, nclusts = nclusts)

    ests = Basu_estimation()
    ests.df= sim.df
    ests.GRM = sim.GRM
    ests.mpheno = ["Y"] 
    
    # Cretae covariates for sites if necessary    
    if nsites > 1:
        cs = ["abcd_site"]
    else :
        cs = None
        
    # Estimate always
    AdjHE = ests.estimate(Method = "AdjHE", npc = [nnpc], covars = cs, mpheno = ["Y"], Naive = False)["h2"][0]
    GCTA_est = ests.estimate(Method = "GCTA", npc = [nnpc], covars = cs, mpheno = ["Y"], Naive = False)["h2"][0]
    
    # Estimate when greater than one site
    if nsites > 1 :
        nAdjHE = ests.estimate(Method = "AdjHE", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = True)["h2"][0]
        AdjHE_RE = ests.estimate(Method = "AdjHE", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = False)["h2"][0]
        try :
            nGCTA = ests.estimate(Method = "GCTA", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = True)["h2"][0]
        except FileNotFoundError :
            print("No estimate since sample size too small for GCTA")
            nGCTA = np.nan

        # hard coded no npcs for the time being should fix!!!
        SWD = ests.estimate(Method = "SWD", npc = [0], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = False)["h2"][0]
        Combat = ests.estimate(Method = "Combat", npc = [0], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = False)["h2"][0]

    else :
        nAdjHE = np.nan
        AdjHE_RE = np.nan
        nGCTA = np.nan
        GCTA = np.nan
        SWD = np.nan
        Combat = np.nan
        
    result = [["GCTA", GCTA_est],
              ["nGCTA", nGCTA],
              ["nAdjHE", nAdjHE],
              ["AdjHE", AdjHE],
              ["SWD", SWD],
              ["Combat", Combat],
              ["AdjHE_RE", AdjHE_RE],
              #["MOM_est",MOM_est]
              ]
    result= pd.DataFrame(result, columns = ["Estimator", "Estimate"])        
        
    # store simulation details
    result["sg"] = sigma[0]
    result["ss"] = sigma[1]
    result["se"] = sigma[2]
    result["nsubjects"] = nsubjects
    result["nsites"] = nsites
    result["theta_alleles"] = theta_alleles
    result["nclusts"] = nclusts
    result["prop_causal"] = prop_causal
    result["nnpc"] = nnpc
    result["nSNPs"] = nSNPs
    result["site_comp"] = site_comp
    result["dominance"] = dominance
    result["site_dep"] = site_dep
    result["site_het"] = site_het
    return result
 
#%%
result = sim_n_est(nsubjects = 1000, sigma = [0.5,0.25, 0.5], site_comp = "EQUAL", nsites = 25,
              theta_alleles =0.9, nclusts =5,  dominance=3, prop_causal=0.25, site_dep=False, nnpc = 1,
              nSNPs=200, phens = 2,  site_het = True)
print(result)

#%%
N = 1000
nc = 2
ns = 2
site_het = False


df = sim_n_est(nsubjects = N, sigma = [0.5,0.25, 0.25], site_comp = "IID", nsites = ns,
              theta_alleles =0.5, nclusts =nc,  dominance=3, prop_causal=0.25, site_dep=False, nnpc = nc,
              nSNPs=20, phens = 2, reps = 10, site_het = False) 
sns.boxplot(data= df, x=  "Estimator", y = "Estimate")
plt.axhline(0.5)
plt.xticks(rotation=45)
plt.title(f"N: {N}, #Clust: {nc}, #Sites: {ns}, Het: {site_het}")

plt.savefig(f"Simulations/Small_verifications/N{N}_C{nc}_S{ns}_Het{site_het}.png")
