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
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
from AdjHE.estimation.all_estimators import Basu_estimation




rng = np.random.default_rng()

def sim_n_est(nsubjects = 1000, sigma = [0.5,0.25, 0.25], site_comp = "IID", nsites = 30,
              theta_alleles =0.5, nclusts =1, dominance=3, prop_causal=0.25, site_dep=False, nnpc = 1,
              nSNPs=20, phens = 2, site_het = False, clusters_differ = False, cov_effect = True,
              ortho_cov = True, random_BS=True) :
    sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
    # Run through full simulation and estimation
    sim.full_sim(nsites= nsites, sigma= sigma, phens = phens, nclusts = nclusts, clusters_differ = clusters_differ,
                 prop_causal = prop_causal, cov_effect = cov_effect, ortho_cov = ortho_cov, random_BS = random_BS)

    ests = Basu_estimation()
    ests.df= sim.df
    ests.GRM = sim.GRM
    ests.mpheno = ["Y1"] 
    
    # Cretae covariates for sites if necessary    
    if cov_effect :
        FE= ["Xc"]
    else :
        FE = []
        
    if nsites >1 :
        FE += ["abcd_site"]
    
    logging.info(sim.df.columns.tolist())

    # List of methods to estimate
    Methods = ["AdjHE", "GCTA"]
    # Meta 
    Metas = [True, False]
    # Treating pcs
    pc_2moments = [True, False]
    random_groups = [None, "abcd_site"]

    full_results = pd.DataFrame({}, columns= ["Method", "Meta", "PC_moment", "RV", "h2", "V(h2)", "time"])
    # iterate over all combinations of Methods, Meta, and pc_2moment 
    for Method, Meta, pc_2nd, random_group in itertools.product(Methods, Metas, pc_2moments, random_groups):
        try :
            r = ests.estimate(Method = Method, npc = [nnpc], fixed_effects= FE, mpheno =ests.mpheno[0],
                        Naive = Meta, pc_2moment = pc_2nd, random_groups = random_group)

            result = pd.DataFrame({"Method" : Method,
                                  "Meta" : Meta,
                                  "PC_moment" : pc_2nd,
                                  "RV" : random_group,
                                  "h2" : r["h2"][0],
                                  "Var(h2)" : r["var(h2)"][0],
                                   "time" : r["time"][0]}, index= [0])
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
            full_results = pd.concat([full_results, result], axis = 0)
            return full_results
        except TypeError :
            return None
        except np.linalg.LinAlgError :
            return None
        
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
