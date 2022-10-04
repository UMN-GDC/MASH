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
from scipy.stats import random_correlation
from functions.load_data import ReadGRMBin, build_grm
from functions.eigenvector_outters import multiple_outer
from functions.simulation_helpers.simulate_GRM_phenos import sim_GRM, simulate_phenotypes

# NUMBER OF REPS FOR EACH CONFIGURATION
reps = 5

#%% create the domain of variance contributions over which to simulate
sg = np.array(range(5)) /3
ss = np.array(range(5))/3
se = np.array(range(5))/3
sigmas = np.array(list((itertools.product(sg,ss, se))))
# only grab the ones that sum to 1 to make interpretation more direct
sigmas = sigmas[sigmas.sum(axis= 1) == 1]

#%% load GRM
G = ReadGRMBin("/panfs/roc/groups/3/rando149/coffm049/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/full/full")
GRM, df = build_grm(G)
del G
n= GRM.shape[0]
#%% Simulate/load a random GRM (nonsparse)

sim_GRM(n, "simulations/Random_corr")

print("loading simulated GRM")
A = np.load("simulations/Random_corr.npy")[0:n, 0:n]
#%% load data on sites

df2= pd.read_csv("/panfs/roc/groups/3/rando149/coffm049/ABCD/Results/02_Phenotypes/Covars.tsv", sep = "\t")

df2 = df.merge(df2, left_on = ["fid", "iid"], right_on = ["FID", "IID"])[["FID","IID", "abcd_site"]]

# Then get the indices to keep for the GRM
GRM_keep = [iid in df2 for iid in df.iid]

GRM = GRM[GRM_keep,:][:, GRM_keep]

#%% Simulate model Y = 0 + e,   e ~ N(0, sgA + ssS + se I)
# with varying variances attached to each covariance structure

df = simulate_phenotypes(GRM, df, sigmas, "EUR_norels_ABCD", reps = 5)
df2 = simulate_phenotypes(A, df, sigmas, "EUR_norels_sim", reps = 5 )


df.to_csv("Simulated_ABCD_full_phenotypes.csv", index= False)
df2.to_csv("Simulated_non_sparse_full_phenotypes.csv", index= False) 
    
