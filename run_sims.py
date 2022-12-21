#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 18:08:00 2022

@author: christian
"""


import os
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from functions.simulation_helpers.Sim_generator import pheno_simulator
from functions.Estimation.all_estimators import Basu_estimation

rng = np.random.default_rng()

#%%

# nsubjects = 100; sigma = [0.5,0.25, 0.25]; site_comp = "IID"
# nsites = 2; theta_alleles =0.5; nclusts =1
# dominance=3; prop_causal=0.25; site_dep=False
# nnpc = 1
# nSNPs=20
# rep = 2
# phens = 2



def sim_n_est(nsubjects = 1000, sigma = [0.5,0.25, 0.25], site_comp = "IID", nsites = 30,
              theta_alleles =0.5, nclusts =1, dominance=3, prop_causal=0.25, site_dep=False, nnpc = 1,
              nSNPs=20, phens = 2, reps = 10, all_ests = True) :
    
    results = pd.DataFrame({})
    for i in range(reps):
        sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
        # Run through full simulation and estimation
        sim.full_sim(nsites= nsites, sigma=sigma, phens = phens)
        
        ests = Basu_estimation()
        ests.df= sim.df
        ests.GRM = sim.GRM
        ests.mpheno = ["Y"]
        if nsites  == 1 :
            cs = None
        else :
            cs = ["abcd_site"]
        ests.looping(covars = cs, npc = [nnpc], mpheno = ["Y"], loop_covars = False)
        
        # Fixed effects AdjHE
        AdjHE_FE = ests.estimate(npc = [nnpc], Method = "AdjHE", Naive = False, covars = True)["h2"][0]
        
        if all_ests :
        
            # SWD
            SWD_est = ests.estimate(npc = [nnpc], Method = "SWD", Naive = False, RV = "abcd_site")["h2"][0]
        
            # COMBAT
            Combat_est = ests.estimate(npc = [nnpc], Method = "Combat", Naive = False, RV = "abcd_site")["h2"][0]
        
            
            # Naive GCTA
            nGCTA_est = ests.estimate(npc = [nnpc], Method = "GCTA", Naive = True, RV = "abcd_site")["h2"][0]
            
            # Naive AdjHE
            nAdjHE_est = ests.estimate(npc = [nnpc], Method = "AdjHE", Naive = True, RV = "abcd_site")["h2"][0]
    
        else :
            SWD_est = np.nan
            Combat_est = np.nan
            nGCTA_est = np.nan
            nAdjHE_est = np.nan
        
        # Standard GCTA
        GCTA_est = ests.estimate(npc = [nnpc], Method = "GCTA", Naive = False)["h2"][0]
        
        # Random effects AdjHE
        ests.looping(covars=  None, npc = [nnpc], mpheno = ["Y"], loop_covars = False)
        AdjHE_RE = ests.estimate(npc = [nnpc], Method = "AdjHE", Naive = False, RV = "abcd_site")["h2"][0]
        
        # MOM estimator
        # MOM_est = ests.estimate(npc = [nnpc], Method = "MOM", Naive = False, covars= ["abcd_site"])["h2"][0]
    
        
        # Make a list of lists with estmiator and estimates
        #result = [[est, eval(est)] for est in [#]]
        result = [["GCTA_est", GCTA_est],
                  ["nGCTA_est", nGCTA_est],
                  ["nAdjHE_est", nAdjHE_est],
                  ["AdjHE_FE", AdjHE_FE],
                  ["SWD_est", SWD_est],
                  ["Combat_est", Combat_est],
                  ["AdjHE_RE", AdjHE_RE],
                  #["MOM_est",MOM_est]
                  ]
        result= pd.DataFrame(result, columns = ["Estimator", "Estimate"])
        # add it to list of results
        results = results.append(result, ignore_index = True)
        
        
    # store simulation details
    results["sg"] = sigma[0]
    results["ss"] = sigma[1]
    results["se"] = sigma[2]
    results["nsubjects"] = nsubjects
    results["nsites"] = nsites
    results["theta_alleles"] = theta_alleles
    results["nclusts"] = nclusts
    results["prop_causal"] = prop_causal
    results["nnpc"] = nnpc
    results["nSNPs"] = nSNPs
    results["site_comp"] = site_comp
    results["dominance"] = dominance
    results["site_dep"] = site_dep
    return results

#%%

def sim_experiment(nsubjectss = [1000], sigmas = [[0.5,0.25, 0.25]], site_comps = ["IID"], nsites = [25],
              theta_alleless = [0.9], nclustss = [5], dominances= [3], prop_causals= [0.05], site_deps= [False], nnpcs = [1],
              nSNPss= [200], phenss= [2], reps = 10, all_ests = True) :
    # Seed empty dataframe
    sim_results = pd.DataFrame()
    
    for (nsubjects, sigma, site_comp, nsite, theta_alleles, nclusts, dominance, prop_causal, site_dep, nnpc, nSNPs, phens
         ) in itertools.product(nsubjectss, sigmas, site_comps, nsites, theta_alleless, nclustss, dominances, prop_causals, site_deps, nnpcs, nSNPss, phenss) :
        
        result = sim_n_est(nsubjects = nsubjects, sigma = sigma, site_comp = site_comp, nsites = nsite,
                           theta_alleles = theta_alleles, nclusts = nclusts, dominance= dominance, prop_causal= prop_causal, 
                           site_dep= site_dep, nnpc = nnpc,
                           nSNPs=nSNPs, phens = phens, reps = reps, all_ests = all_ests)
        sim_results= sim_results.append(result, ignore_index = True)
    
    return sim_results
                                                     


#%% Create domain for simulations
# sgs = [0, 0.25, 0.5]
# sss = [0, 0.25]
# ses = [0, 0.25, 0.5, 0.75, 1]

# sgs = [0.5]
# sss = [0.25]
# ses = [0.25]

# sigmas = []
# for sg, ss, se in itertools.product(sgs, sss, ses) :
#     if sg + ss + se == 1 :
#         if sg != 0 :
#             sigmas += [[sg, ss, se]]
#         elif (sg ==0) and (ss == 0) :
#             sigmas += [[sg, ss, se]]
#%%  
# N  = 1000
# ns = 25
# nc = 5
# # #%%
# df = sim_experiment(nsubjectss= [N], reps= 25, nsites=[ns], site_comps = ["EQUAL"], sigmas = sigmas, nnpcs = [5], nclustss=[nc],
#                     all_ests = False)
# df.to_csv("Simulations/Sim_working_Combat1.csv", header=  True, index= False)

# g = sns.FacetGrid(df, col="sg",  row="ss", sharey = False)
# g.map(sns.boxplot, "Estimator", "Estimate")
# 
#%%
# g = sns.boxplot(data = df, x= "Estimator", y="Estimate")
# g.axhline(0.5)
# plt.xticks(rotation=45)
# plt.title("N = " + str(N) + ", #sites= " + str(ns) + ", #clusts=" + str(nc))

# df[["Estimator", "Estimate"]].groupby("Estimator").agg(['mean', 'std'])
