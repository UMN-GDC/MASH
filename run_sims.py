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

def sim_experiment(nsubjectss = [1000], sigmas = [[0.5,0.25, 0.25]], site_comps = ["IID"], nsites = [25],
              theta_alleless = [0.9], nclustss = [5], dominances= [3], prop_causals= [0.05], site_deps= [False], nnpcs = [1],
              nSNPss= [200], phenss= [2], reps = 10, site_het = False) :
    # Seed empty dataframe
    sim_results = pd.DataFrame()
    
    for (nsubjects, sigma, site_comp, nsite,
         theta_alleles, nclusts, dominance, prop_causal,
         site_dep, nnpc, nSNPs, phens,
         __
         ) in itertools.product(nsubjectss, sigmas, site_comps, nsites,
                                theta_alleless, nclustss, dominances, prop_causals,
                                site_deps, nnpcs, nSNPss, phenss,
                                range(reps)) :
                                
        nnpc = 2 * nclusts
                                
        result = sim_n_est(nsubjects = nsubjects, sigma = sigma, site_comp = site_comp, nsites = nsite,
                           theta_alleles = theta_alleles, nclusts = nclusts, dominance= dominance, prop_causal= prop_causal, 
                           site_dep= site_dep, nnpc = nnpc,
                           nSNPs=nSNPs, phens = phens, site_het = site_het)
        sim_results= sim_results.append(result, ignore_index = True)
        # Remove any temps
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
