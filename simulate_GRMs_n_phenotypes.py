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
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
import itertools
import numpy as np
import pandas as pd
from scipy.stats import random_correlation
from functions.load_data import ReadGRMBin, build_grm
from functions.eigenvector_outters import multiple_outer
from functions.simulation_helpers.simulate_GRM_phenos import sim_GRM, simulate_phenotypes

#%% create the domain of variance contributions over which to simulate
sg = np.array(range(5)) /5
ss = np.array(range(5))/5
se = np.array(range(5))/5
sigmas = np.array(list((itertools.product(sg,ss, se))))
# only grab the ones that sum to 1 to make interpretation more direct
sigmas = sigmas[sigmas.sum(axis= 1) == 1]

#%% load GRM
G = ReadGRMBin("simulations/Data/EUR_no_rels")
GRM, df = build_grm(G)
del G
n= GRM.shape[0]
#%% Simulate/load a random GRM (nonsparse)

sim_GRM(11878, "simulations/Random_corr.npy")

A = np.load("simulations/Random_corr.npy")[0:n, 0:n]
#%% load data on sites

df2= pd.read_csv("simulations/Data/Covars.tsv", sep = " ")

df = df.merge(df2, left_on = ["fid", "iid"], right_on = ["FID", "IID"])[["FID","IID", "abcd_site"]]



#%% Simulate model Y = 0 + e,   e ~ N(0, sgA + ssS + se I)
# with varying variances attached to each covariance structure

df = simulate_phenotypes(GRM, df, sigmas, "EUR_no_rels")

    
