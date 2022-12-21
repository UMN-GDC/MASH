#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 16:23:04 2022

@author: christian
"""
from functions.simulation_helpers.Sim_generator import pheno_simulator
from run_sims import sim_experiment, Basu_estimation
import matplotlib.pyplot as plt
import seaborn as sns

#%%

# sim = pheno_simulator(nsubjects= 2000, nSNPs = 200)
# # Run through full simulation and estimation
# sim.full_sim(nsites= 2, sigma=[ 0.5, 0.01, 0.5], phens = 1, nclusts = 2)

# ests = Basu_estimation()
# ests.df= sim.df
# ests.GRM = sim.GRM
# ests.mpheno = ["Y"]
# ests.looping(covars = None, npc = [2], mpheno = ["Y"], loop_covars = False)

# ests.estimate(Method = "AdjHE", npc = [1], covars= None)
#%%
N = 1000
nc = 3
ns = 2
site_het = True


df = sim_experiment(nsubjectss = [N], sigmas = [[0.5,0.25, 0.25]], site_comps = ["EQUAL"], nsites = [ns],
              theta_alleless = [0.9], nclustss = [nc], dominances= [3], prop_causals= [0.05], site_deps= [False], nnpcs = [nc],
              nSNPss= [200], phenss= [2], reps = 10, all_ests = False, site_het = site_het)


sns.boxplot(data= df, x=  "Estimator", y = "Estimate")
plt.axhline(0.5)
plt.xticks(rotation=45)
plt.title(f"N: {N}, #Clust: {nc}, #Sites: {ns}, Het: {site_het}")

plt.savefig(f"Simulations/Small_verifications/N{N}_C{nc}_S{ns}_Het{site_het}.png")

