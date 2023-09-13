#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 18:08:00 2022

@author: christian
"""

import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from tqdm.auto import tqdm
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
from Estimate.estimators.all_estimators import h2Estimation




rng = np.random.default_rng()

def sim_n_est(nsubjects = 1000, h2 = 0.5, nsites = 30,
              nclusts =1, nnpc = 1,
              nSNPs=20, phens = 2,
              random_BS=True) :
    sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
    # Run through full simulation and estimation
    sim.full_sim(nsites= nsites, h2= h2, phens = phens, nclusts = nclusts, random_BS = random_BS)

    ests = h2Estimation()
    ests.df= sim.df
    ests.GRM = sim.GRM
    ests.mpheno = ["Y1"] 
    
    FE= ["Xc"]
        
    if nsites >1 :
        FE += ["abcd_site"]
    
    logging.info(sim.df.columns.tolist())

    full_results = pd.DataFrame({}, columns= ["Method", "Meta", "RV", "h2", "V(h2)", "time"])
    # iterate over all combinations of Methods and Meta 
    for (Met, M, rg) in itertools.product(["AdjHE", "GCTA"], [True, False], [None, "abcd_site"]):
        
        try :
            r = ests.estimate(Method = Met, npc = [nnpc], fixed_effects= FE, mpheno =ests.mpheno[0],
                        Naive = M, random_groups = rg)
            result = pd.DataFrame({"Method" : Met,
                                  "Meta" : M,
                                  "RV" : rg,
                                  "h2" : r["h2"][0],
                                  "Var(h2)" : r["var(h2)"][0],
                                   "time" : r["time"][0],
                                  "sg" : sigma[0],
                                  "ss" : sigma[1],
                                  "se" : sigma[2],
                                  "nsubjects" : nsubjects,
                                  "nsites" : nsites,
                                  "theta_alleles" : theta_alleles,
                                   "nclusts" : nclusts,
                                   "prop_causal" : prop_causal,
                                   "nnpc" : nnpc,
                                   "nSNPs" : nSNPs, 
                                   "site_comp" : site_comp,
                                   "dominance" : dominance,
                                   "site_dep" : site_dep,
                                   "site_het" : site_het}, index= [0])
            full_results = pd.concat([full_results, result], ignore_index= True)
        except TypeError :
            pass
        except np.linalg.LinAlgError :
            pass

        return full_results
        
def sim_experiment(nsubjectss = [1000], sigmas = [[0.5,0.25, 0.25]], site_comps = ["IID"], nsites = [25],
              theta_alleless = [0.9], nclustss = [5], dominances= [3], prop_causals= [0.05], site_deps= [False], nnpcs = [1],
              nSNPss= [200], phenss= [2], reps = 10, site_het = False, clusters_differ = False, cov_effect = True,
              ortho_cov = True, random_BS = True) :
    # Seed empty dataframe
    sim_results = pd.DataFrame()
    
    for (nsubjects, sigma, site_comp, nsite,
         theta_alleles, nclusts, dominance, prop_causal,
         site_dep, nnpc, nSNPs, phens,
         __
         ) in tqdm(itertools.product(nsubjectss, sigmas, site_comps, nsites,
                                theta_alleless, nclustss, dominances, prop_causals,
                                site_deps, nnpcs, nSNPss, phenss,
                                range(reps)), desc = "Simulation progress") :
                                
        # nnpc = 2 * nclusts
                                
        result = sim_n_est(nsubjects = nsubjects, sigma = sigma, site_comp = site_comp, nsites = nsite,
                           theta_alleles = theta_alleles, nclusts = nclusts, dominance= dominance, prop_causal= prop_causal, 
                           site_dep= site_dep, nnpc = nnpc,
                           nSNPs=nSNPs, phens = phens, site_het = site_het, clusters_differ = clusters_differ, cov_effect = cov_effect,
                           ortho_cov = ortho_cov, random_BS = random_BS)
        sim_results = pd.concat([sim_results, result], ignore_index=True)
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
#N  = 500
#ns = 5
#nc = 2
#sigmas = [[0.5,0.25,0.25]]
#%%

#df = sim_experiment(nsubjectss= [N], reps= 5, nsites=[1], site_comps = ["EQUAL"], sigmas = sigmas, nnpcs = [1], nclustss=[nc], 
#                    ortho_cov = True, cov_effect= True, phenss= [3], random_BS = False)
#logging.info(df[["GCTA", "AdjHE", "AdjHE_RE", "Combat", "SWD"]].mean())

#%%
# # df.to_csv("Simulations/Sim_working_Combat1.csv", header=  True, index= False)

# g = sns.FacetGrid(df, col="sg",  row="ss", sharey = False)
# g.map(sns.boxplot, "Estimator", "Estimate")
# # 

#g = sns.boxplot(data = df, x= "Estimator", y="Estimate")
# g.axhline(0.5)
# plt.xticks(rotation=45)
# plt.title("N = " + str(N) + ", #sites= " + str(ns) + ", #clusts=" + str(nc))

# df[["Estimator", "Estimate"]].groupby("Estimator").agg(['mean', 'std'])
