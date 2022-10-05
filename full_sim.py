#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:32:40 2022

@author: christian
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate a GRM (nonsparse) and associated phenotypes
Simulate phentotypes off of the ABCD GRM
All of these will be simulated with varying degrees of contributino from the GRM, sites, and noise
Created on Fri Sep 23 10:00:54 2022
s
@author: christian
"""
import os 
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
os.chdir("/panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit")
import itertools
import numpy as np
import pandas as pd
from functions.load_data import ReadGRMBin, build_grm
from functions.simulation_helpers.simulate_GRM_phenos import sim_GRM, simulate_phenotypes

# NUMBER OF REPS (reps) FOR EACH CONFIGURATION and number of subejcts (n)
reps = 3
steps = 3 
n = 5000

#%% load GRM
GRM, ids = build_grm(ReadGRMBin("/panfs/roc/groups/3/rando149/coffm049/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/full/full"))
GRM = GRM[0:n,:][:, 0:n]
ids = ids.iloc[0:n,:]

#%% load covariates 
df = pd.read_csv("/panfs/roc/groups/3/rando149/coffm049/ABCD/Results/02_Phenotypes/Covars.tsv", sep = "\t")
# use the order from the ids data to match the order of the GRM
df = ids.merge(df, left_on = ["fid", "iid"], right_on = ["FID", "IID"])[["FID","IID", "abcd_site"]]

# Then get the indices to keep for the GRM
GRM_keep = [i in list(df.IID) for i in ids.iid]
GRM = GRM[GRM_keep,:][:, GRM_keep]

n = GRM.shape[0]

#%% Simulate/load a random GRM (nonsparse)
sim_GRM(n, "simulations/Random_corr.npy")
print("loading simulated GRM")
A = np.load("simulations/Random_corr.npy")[0:n, 0:n]

#%% Simulate model Y = 0 + e,   e ~ N(0, sgA + ssS + se I)
# with varying variances attached to each covariance structure

#%% create the domain of variance contributions over which to simulate
sg = np.array(range(5)) /steps
ss = np.array(range(5))/steps
se = np.array(range(5))/steps
sigmas = np.array(list((itertools.product(sg,ss, se))))
# only grab the ones that sum to 1 to make interpretation more direct
sigmas = sigmas[sigmas.sum(axis= 1) == 1]

# And simulated values to the datframe then save it
df = simulate_phenotypes(GRM, df, sigmas, "ABCD", reps = reps)
df = simulate_phenotypes(A, df, sigmas, "synth", reps = reps )
df.to_csv("Simulated_full_phenotypes.csv", index= False, sep = "\t")
    
