#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scripts for making figures of the estimation process using this tool

Created on Thu Dec 15 13:01:26 2022

@author: christian
"""

import os 
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
from functions.simulation_helpers.Sim_generator import pheno_simulator
from functions.Estimation.all_estimators import Basu_estimation
import seaborn as sns

#%%

sim = pheno_simulator(nsubjects = 5000, nSNPs = 5000)
sim.full_sim(sigma = [0.5, 0.25,0.25], nclusts = 5, nsites = 30, theta_alleles=0.9)


ests = Basu_estimation()
ests.df= sim.df
ests.GRM = sim.GRM

#%%

ests.GRM_vis(sort_by = "subj_ancestries")
ests.GRM_vis(sort_by = "subj_ancestries",   npc=4)
#%%
trans = PCA(n_components=2).fit_transform(ests.GRM)
#%%
ax = sns.scatterplot(x=trans[:, 0], y=trans[:, 1],
                hue=ests.df.subj_ancestries.astype(str),
                hue_order= [str(i) for i in range(len(set(ests.df.subj_ancestries)))])
ax.set(xlabel='PC1', ylabel='PC2')



