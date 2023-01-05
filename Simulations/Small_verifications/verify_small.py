#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 16:23:04 2022

@author: christian
"""
import os
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
from functions.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pandas as pd
from functions.Estimation.all_estimators import Basu_estimation
from run_sims import  sim_n_est, sim_experiment
import matplotlib.pyplot as plt
import seaborn as sns
import itertools

#%%
N = 1000
nc = 5
ns = 5
site_het = False


df = sim_n_est(nsubjects = N, sigma = [0.5,0.01, 0.5], site_comp = "IID", nsites = ns,
              theta_alleles =0.5, nclusts =nc,  dominance=3, prop_causal=0.25, site_dep=False, nnpc = 0,
              nSNPs=20, phens = 2, site_het = False) 
print(df)

#%%
sgs = [0, 0.25, 0.5]
sss = [0, 0.25]
ses = [0, 0.25, 0.5, 0.75, 1]

sgs = [0.5]
sss = [0.25]
ses = [0.25]

sigmas = []
for sg, ss, se in itertools.product(sgs, sss, ses) :
    if sg + ss + se == 1 :
        if sg != 0 :
            sigmas += [[sg, ss, se]]
        elif (sg ==0) and (ss == 0) :
            sigmas += [[sg, ss, se]]


N  = 1000
ns = 25
nc = 5
site_het = True
# #%%
df = sim_experiment(nsubjectss= [N], reps= 10, nsites=[ns], site_comps = ["EQUAL"], sigmas = sigmas, nnpcs = [nc], nclustss=[nc],
                    site_het = site_het)

sns.boxplot(data= df, x=  "Estimator", y = "Estimate")
plt.axhline(0.66)
plt.xticks(rotation=45)
plt.title(f"N: {N}, #Clust: {nc}, #Sites: {ns}, Het: {site_het}")

plt.savefig(f"Simulations/Small_verifications/N{N}_C{nc}_S{ns}_Het{site_het}.png")


#%%
ns = 5
nc = 4
N = 500
site_het = False

sim= pheno_simulator(nsubjects=N, nSNPs=100)
sim.full_sim(sigma = [0.5,0.25,0.25], nsites = ns, nclusts = nc, site_het = site_het)


ests = Basu_estimation()
ests.df= sim.df
ests.GRM = sim.GRM
ests.mpheno = ["Y"] 

# Cretae covariates for sites if necessary    
if ns > 1:
    cs = ["abcd_site"]
else :
    cs = None


GCTA_est = ests.estimate(Method = "AdjHE", npc = [0], covars = cs, mpheno = ["Y"], RV = None, Naive = False)
print(GCTA_est)

#%%


AdjHE = ests.estimate(Method = "AdjHE", npc = [nc -1 ], covars = cs, mpheno = ["Y"], Naive = False)["h2"][0]
print(AdjHE)
#%%

ests2 = Basu_estimation()
ests2.df = sim.df.loc[sim.df["abcd_site"] == 1,:]
ests2.GRM = sim.GRM[sim.df["abcd_site"] == 1,:][:,sim.df["abcd_site"] == 1]
ests.mpheno= ["Y"]
if ns > 1:
    cs = ["abcd_site"]
else :
    cs = None

GCTA_est= ests2.estimate(Method = "GCTA", npc = [0 ], covars = [], mpheno = ["Y"], RV = None, Naive = False)
print(GCTA_est)
