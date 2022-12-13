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
import itertools
from functions.simulation_helpers.Sim_generator import pheno_simulator
from functions.Estimation.all_estimators import Basu_estimation

rng = np.random.default_rng(123)

#%%

nsubjects = 100; sigma = [0.5,0.25, 0.25]; site_comp = "IID"
nsites = 2; theta_alleles =0.5; nclusts =1
dominance=3; prop_causal=0.25; site_dep=False
nnpc = 1
nSNPs=20
rep = 2
phens = 2


def sim_experiment(nsubjectss = [100], site_comps=["IID"], nSNPss = [20],
                    nsitess=[30], theta_alleless=[0.5], nclustss=[5], dominances=[5],
                    prop_causals=[0.7], site_deps=[False], reps=25,
                    nnpcs=[0], sigmas = [[0.5,0.25,0.25]], phenss = [2]) :
    # Seed empty dataframe
    sim_results = pd.DataFrame()
    # loop over all simulation variables and number of repetitions
    for nsubjects, sigma, site_comp, nsites, theta_alleles, nclusts, dominance, prop_causal, site_dep, nnpc, __, nSNPs, phens in itertools.product(nsubjectss, sigmas, site_comps, nsitess,
                                                                                                                theta_alleless, nclustss, dominances, 
                                                                                                                prop_causals, site_deps, 
                                                                                                                nnpcs, range(reps), nSNPss, phenss) :
        sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
        # Run through full simulation and estimation
        sim.full_sim(nsites= nsites, sigma=sigma, phens = phens)
        
        ests = Basu_estimation()
        ests.df= sim.df
        ests.GRM = sim.GRM
        ests.mpheno = ["Y"]
        ests.looping(covars = ["abcd_site"], npc = [nnpc], mpheno = ["Y"], loop_covars = False)
        
        
        # Standard GCTA
        GCTA_est = ests.estimate(npc = [nnpc], Method = "GCTA", Naive = False)["h2"][0]
        
        # Naive GCTA
        nGCTA_est = ests.estimate(npc = [nnpc], Method = "GCTA", Naive = True, RV = "abcd_site")["h2"][0]
        
        # Naive AdjHE
        nAdjHE_est = ests.estimate(npc = [nnpc], Method = "AdjHE", Naive = True, RV = "abcd_site")["h2"][0]
        # Fixed effects AdjHE
        AdjHE_FE = ests.estimate(npc = [nnpc], Method = "AdjHE", Naive = False, covars = True)["h2"][0]
        # SWD
        SWD_est = ests.estimate(npc = [nnpc], Method = "SWD", Naive = False, RV = "abcd_site")["h2"][0]

        # COMBAT
        Combat_est = ests.estimate(npc = [nnpc], Method = "Combat", Naive = False, RV = "abcd_site")["h2"][0]

        
        # Random effects AdjHE
        ests.looping(covars=  None, npc = [nnpc], mpheno = ["Y"], loop_covars = False)
        AdjHE_RE = ests.estimate(npc = [nnpc], Method = "AdjHE", Naive = False, RV = "abcd_site")["h2"][0]
        
        # MOM estimator
        # MOM_est = ests.estimate(npc = [nnpc], Method = "MOM", Naive = False, covars= ["abcd_site"])["h2"][0]

        
        # Make a list of lists with estmiator and estimates
        result = [[est, eval(est)] for est in ["GCTA_est", "nGCTA_est", "nAdjHE_est", 
                                               "AdjHE_FE", "SWD_est", "Combat_est", 
                                               "AdjHE_RE"#, "MOM_est"
                                               ]]
       
        
        # Make it a dataframe
        result = pd.DataFrame(result, columns = ["Estimator", "Estimate"])
        
        # add simulation details
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
        
        sim_results = pd.concat([sim_results, result], ignore_index=True)
        
    # Return results in a tall dataframe
    return sim_results

#%%


#def plot_save(results, out) :
#    results["h2"] = results.sg / (results.sg + results.ss + results.se)
#    results["h2"][np.isnan(results.h2)] = 0
#    results.to_csv(out + ".csv")

    # Plot results and store at specified locaation
    # fig = px.violin(results, x="variable", y="value", color="variable", facet_col="h2", facet_row = "ss")
#    fig = px.box(results, x="variable", y="value",
#                 color="variable", facet_col="h2", facet_row="ss")
#    fig.update_yaxes(matches=None)
#    fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
#    fig.update_xaxes(matches=None)
#    fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True))
#    fig.update_layout(
#        font=dict(size=20)
#    )
#    plot(fig, filename=out + ".html")

#%% Create domain for simulations
sgs = [0, 0.25, 0.5]
sss = [0, 0.25]
ses = [0, 0.25, 0.5, 0.75, 1]

sigmas = []
for sg, ss, se in itertools.product(sgs, sss, ses) :
    if sg + ss + se == 1 :
        if sg != 0 :
            sigmas += [[sg, ss, se]]
        elif (sg ==0) and (ss == 0) :
            sigmas += [[sg, ss, se]]
            
#%%
cool_results = sim_experiment(nsubjectss = [100], site_comps=["IID"], nSNPss = [20],
                    nsitess=[2], theta_alleless=[0.75], nclustss=[5], dominances=[5],
                    prop_causals=[0.25], site_deps=[False], reps=2,
                    nnpcs=[5], sigmas = sigmas, phenss = [2])

#%%

cool_results.to_csv("Simulations/Sim_results2.csv", header=  True, index= False)

