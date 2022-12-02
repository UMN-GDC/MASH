#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 18:08:00 2022

@author: christian
"""


import os
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")

import numpy as np
import pandas as pd
import itertools
from functions.simulation_helpers.Sim_generator import pheno_simulator
from functions.Estimation.all_estimators import Basu_estimation


rng = np.random.default_rng(123)

#%%


def sim_experiment(nsubjectss = [100], site_comps=["IID"], nSNPss = [20],
                    nsitess=[30], theta_alleless=[0.5], nclustss=[5], dominances=[5],
                    prop_causals=[0.7], site_deps=[False], reps=25,
                    nnpcs=[0], sigmas = [[0.5,0.25,0.25]]) :
    # Seed empty dataframe
    sim_results = pd.DataFrame()
    # loop over all simulation variables and number of repetitions
    for nsubjects, sigma, site_comp, nsites, theta_alleles, nclusts, dominance, prop_causal, site_dep, nnpc, __, nSNPs in itertools.product(nsubjectss, sigmas, site_comps, nsitess,
                                                                                                                theta_alleless, nclustss, dominances, 
                                                                                                                prop_causals, site_deps, 
                                                                                                                nnpcs, range(reps), nSNPss) :
        sim = pheno_simulator(nsubjects= nsubjects, nSNPs = nSNPs)
        # Run through full simulation and estimation
        sim.full_sim(nsites= nsites, sigma= sigma, )
        
        ests = Basu_estimation()
        ests.df= sim.df
        ests.GRM = sim.GRM
        ests.mpheno = ["Y"]
        ests.looping(covars = ["abcd_site"], npc = nnpcs, mpheno = ["Y"], loop_covars = False)
        
    # specify parameters and estimate  column names
    params = ['sg', 'ss', 'se', 'nsub','nsites', 'nclusts', 'prop_causal', 'site_comp',
              'nSNPs', 'theta_alleles', "dominance", 'npc']
    ests = ['AdjHE_FE', 'AdjHE_RE', 'Naive_AdjHE', 'GCTA', 'Naive_GCTA']
    # Make results tall
    sim_results = (pd.DataFrame(sim_results)
                .melt(id_vars= params, value_vars= ests)
                )
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



sim = pheno_simulator(nsubjects= 1000, nSNPs = 50)
# Run through full simulation and estimation
sim.full_sim(nsites= 30, sigma=[0.5, 0.25, 0.25], )

ests = Basu_estimation()
ests.df= sim.df
ests.GRM = sim.GRM
ests.mpheno = ["Y"]
ests.looping(covars = ["abcd_site"], npc = [5], mpheno = ["Y"], loop_covars = False)
#%%

# Standard GCTA
ests.estimate(npc = [5], Method = "GCTA", Naive = False)
# Naive GCTA
ests.estimate(npc = [5], Method = "GCTA", Naive = True, RV = "abcd_site")

# Naive AdjHE
ests.estimate(npc = [5], Method = "AdjHE", Naive = True, RV = "abcd_site")
# Fixed effects AdjHE
tt = ests.estimate(npc = [5], Method = "AdjHE", Naive = False, covars = True)
# Random effects AdjHE
ests.estimate(npc = [5], Method = "AdjHE", Naive = False, RV = "abcd_site")



