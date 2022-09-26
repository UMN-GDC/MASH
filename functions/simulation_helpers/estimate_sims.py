#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make estimations on both types of simulations ABCD GRM and simulated GRM
Created on Fri Sep 23 12:10:25 2022

@author: christian
"""

import os
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
import numpy as np
import pandas as pd
from functions.load_data import ReadGRMBin, build_grm
from functions.eigenvector_outters import multiple_outer
from functions.two_step_M_est import two_level_regression
import seaborn as sns

#%% load GRM
G = ReadGRMBin("simulations/Data/EUR_no_rels")
GRM, df = build_grm(G)
del G
n = GRM.shape[0]
# load in as many observations as for the ABCD GRM to make it comparable
A = np.load("simulations/Random_corr.npy")[0:n, 0:n]
# Assign each subject to a site
df = pd.read_csv("simulations/just_REs.csv")
df["Zs"] = pd.Series(np.repeat(list(range(10)), 500)[0:n])
df.Zs=  df.Zs.astype(str)


#%% Create the S matrix
Zs = pd.get_dummies(Zs)
Zs = multiple_outer(Zs, Zs)
#%%

h1 = two_level_regression(data= df, GRM = A, fixed = "Zs", dep_var = "sim_820", RV = None)
h2 = two_level_regression(data= df, GRM = GRM, fixed = "Zs", dep_var = "abcd_sim_820", RV = None)
print(h1.params)
print(h2.params)
#%%
h2s = {}
for c in df.columns:
    if "abcd" in c :
        h = two_level_regression(data= df, GRM = GRM, fixed = "Zs", dep_var = c, RV = None).params["GRM"]
        h2s[c] = h
    elif "sim" in c: 
        h = two_level_regression(data= df, GRM = A, fixed = "Zs", dep_var = c, RV = None).params["GRM"]
        h2s[c] = h
    



#%%

ests = pd.DataFrame({"Simulation" : h2s.keys(),
                     "Estimates" : h2s.values()})
ests=  ests.assign(ABCD_GRM = ests.Simulation.str.contains("abcd"),
            h2 = round(ests.Simulation.str.extract('(\d+)').astype(int), -2) / 1000)

#%%
sns.scatterplot(x = "h2", y = "Estimates", hue = "ABCD_GRM", data = ests, s = 10)
plt.plot([0,1], [0,1], lw =0.5)