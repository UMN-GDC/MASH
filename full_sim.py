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
# os.chdir("/panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit")
import itertools
import numpy as np
import pandas as pd
from functions.load_data import ReadGRMBin
from functions.simulation_helpers.simulate_GRM_phenos import sim_GRM, simulate_phenotypes

# NUMBER OF REPS (reps) FOR EACH CONFIGURATION and number of subejcts (n)
reps = 10 
n = 5000

#%% load GRM
prefix = "/panfs/roc/groups/3/rando149/coffm049/ABCD/Results/01_Gene_QC/filters/filter1/GRMs/full/full"
GRM = ReadGRMBin(prefix)
GRM = GRM[0:n,:][:, 0:n]
ids = ids = pd.DataFrame(np.loadtxt(prefix + ".grm.id", delimiter = '\t', dtype = str), columns = ["fid", "iid"])
ids = ids.iloc[0:n,:]

#%% load covariates 
df = pd.read_csv("/panfs/roc/groups/3/rando149/coffm049/ABCD/Results/02_Phenotypes/Covars.tsv", sep = " ")
# use the order from the ids data to match the order of the GRM
df = ids.merge(df, left_on = "iid", right_on = "IID")[["FID","IID", "abcd_site"]]

# Then get the indices to keep for the GRM
GRM_keep = [i in list(df.IID) for i in ids.iid]
GRM = GRM[GRM_keep,:][:, GRM_keep]

n = GRM.shape[0]

#%% Simulate/load a random GRM (nonsparse)
sim_GRM(n, "simulations/phenotypes/synth")
print("loading simulated GRM")
A = ReadGRMBin("simulations/phenotypes/synth")[0:n,0:n]
# Standardize A
A = (A - A.mean(axis = 1))/ A.std(axis = 1)
#%% Simulate model Y = 0 + e,   e ~ N(0, sgA + ssS + se I)
# with varying variances attached to each covariance structure

#%% create the domain of variance contributions over which to simulate
sg = np.array([0, 0.3, 0.6] )
ss = np.array([0, 0.2])
gints= np.array([0, 0.2])
sigmas = np.array(list((itertools.product(sg, ss, gints))))
# add error such that th variances add up to one
sigmas= np.insert(sigmas, 3, 1- np.sum(sigmas,axis = 1), axis=1)
# grab interesting combinations
sigmas = sigmas[[0,2,4, 5, 6, 7, 8, 9, 10, 11] ,:]


# And simulated values to the datframe then save it
df1 = simulate_phenotypes(GRM, df, sigmas, "ABCD", reps = reps)
df2 = simulate_phenotypes(A, df, sigmas, "synth", reps = reps )
df1.drop("abcd_site", axis = 1).to_csv("simulations/phenotypes/ABCD.phen", index= False, sep = " ")
df2.drop("abcd_site", axis = 1).to_csv("simulations/phenotypes/Synth.phen", index= False, sep = " ")

    
