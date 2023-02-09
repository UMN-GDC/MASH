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
from tqdm.auto import tqdm
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
              nSNPs=20, phens = 2, site_het = False, races_differ = False, cov_effect = True,
              ortho_cov = True) :
    sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
    # Run through full simulation and estimation
    sim.full_sim(nsites= nsites, sigma= sigma, phens = phens, nclusts = nclusts, races_differ = races_differ,
                 prop_causal = prop_causal, cov_effect = cov_effect, ortho_cov = ortho_cov)

    ests = Basu_estimation()
    ests.df= sim.df
    ests.GRM = sim.GRM
    ests.mpheno = ["Y"] 
    
    # Cretae covariates for sites if necessary    
    if cov_effect :
        cs = ["Xc"]
    if nsites > 1:
        cs += ["abcd_site"]

    try :
        cs 
    except NameError : 
        cs = []
        
    # Estimate always
    AdjHE = ests.estimate(Method = "AdjHE", npc = [nnpc], covars = cs, mpheno = ["Y"], Naive = False)
    GCTA_est = ests.estimate(Method = "GCTA", npc = [nnpc], covars = cs, mpheno = ["Y"], Naive = False)
    
    # Estimate when greater than one site
    if nsites > 1 :
        nAdjHE = ests.estimate(Method = "AdjHE", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = True)
        AdjHE_RE = ests.estimate(Method = "AdjHE", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = False)
        try :
            nGCTA = ests.estimate(Method = "GCTA", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = True)
        except FileNotFoundError :
            print("No estimate since sample size too small for GCTA")
            nGCTA = np.nan

        SWD = ests.estimate(Method = "SWD", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = False)
        Combat = ests.estimate(Method = "Combat", npc = [nnpc], covars = cs, mpheno = ["Y"], RV = "abcd_site", Naive = False)

    else :
        nAdjHE = pd.DataFrame({"h2" : np.nan, "var(h2)" : np.nan, "Analysis time" : np.nan}, index = [0])
        AdjHE_RE = pd.DataFrame({"h2" : np.nan, "var(h2)" : np.nan, "Analysis time" : np.nan}, index = [0])
        nGCTA = pd.DataFrame({"h2" : np.nan, "var(h2)" : np.nan, "Analysis time" : np.nan}, index = [0])
        SWD = pd.DataFrame({"h2" : np.nan, "var(h2)" : np.nan, "Analysis time" : np.nan}, index = [0])
        Combat = pd.DataFrame({"h2" : np.nan, "var(h2)" : np.nan, "Analysis time" : np.nan}, index = [0])
        
    result = {"GCTA" : GCTA_est["h2"][0],
              "var_GCTA" : GCTA_est["var(h2)"][0],
              "time_GCTA" : GCTA_est["Analysis time"][0],
              
              "nGCTA": nGCTA["h2"][0],
              "var_nGCTA" : nGCTA["var(h2)"][0],
              "time_nGCTA" : nGCTA["Analysis time"][0],

              "nAdjHE": nAdjHE["h2"][0],
              "var_nAdjHE" : nAdjHE["var(h2)"][0],
              "time_nAdjHE" : nAdjHE["Analysis time"][0],

              "AdjHE": AdjHE["h2"][0],
              "var_AdjHE" : AdjHE["var(h2)"][0],
              "time_AdjHE" : AdjHE["Analysis time"][0],

              "SWD": SWD["h2"][0],
              "var_SWD" : SWD["var(h2)"][0],
              "time_SWD" : SWD["Analysis time"][0],
              
              "Combat": Combat["h2"][0],
              "var_Combat" : Combat["var(h2)"][0],
              "time_Combat" : Combat["Analysis time"][0],

              "AdjHE_RE" : AdjHE_RE["h2"][0],
              "var_AdjHE_RE" : AdjHE_RE["var(h2)"][0],
              "time_AdjHE_RE" : AdjHE_RE["Analysis time"][0],
              #["MOM_est",MOM_est]
              }
        
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
    
    result= pd.DataFrame(result, index = [0])        

    return result

#%%

def sim_experiment(nsubjectss = [1000], sigmas = [[0.5,0.25, 0.25]], site_comps = ["IID"], nsites = [25],
              theta_alleless = [0.9], nclustss = [5], dominances= [3], prop_causals= [0.05], site_deps= [False], nnpcs = [1],
              nSNPss= [200], phenss= [2], reps = 10, site_het = False,races_differ = False, cov_effect = True,
              ortho_cov = True) :
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
                           nSNPs=nSNPs, phens = phens, site_het = site_het, races_differ = races_differ, cov_effect = cov_effect,
                           ortho_cov = ortho_cov)
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
N  = 1000
ns = 25
nc = 1
sigmas = [[0.5,0.25,0.25]]
#%%
df = sim_experiment(nsubjectss= [N], reps= 5, nsites=[ns], site_comps = ["EQUAL"], sigmas = sigmas, nnpcs = [nc], nclustss=[nc], ortho_cov = False)
# # df.to_csv("Simulations/Sim_working_Combat1.csv", header=  True, index= False)

# g = sns.FacetGrid(df, col="sg",  row="ss", sharey = False)
# g.map(sns.boxplot, "Estimator", "Estimate")
# # 

#g = sns.boxplot(data = df, x= "Estimator", y="Estimate")
# g.axhline(0.5)
# plt.xticks(rotation=45)
# plt.title("N = " + str(N) + ", #sites= " + str(ns) + ", #clusts=" + str(nc))

# df[["Estimator", "Estimate"]].groupby("Estimator").agg(['mean', 'std'])
