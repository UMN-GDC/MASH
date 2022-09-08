#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 12:13:07 2022

@author: christian
"""

import numpy as np
import pandas as pd
from functions.load_data import load_everything
from functions.parser import read_flags

c_args= {}
c_args['argfile'] = "Example/Argfile.json"
args = read_flags(c_args)

#%% Read in all data
(df, covariates, phenotypes, A, ids) = load_everything(prefix = args["prefix"],
                                                                    pheno_file = args["pheno"], 
                                                                    cov_file= args["covar"], 
                                                                    PC_file= args["PC"],
                                                                    k= args["k"])
#%%

df = pd.read_csv("/home/christian/Research/Stat_gen/tools/Basu_herit/simulations/sim_data.csv")


#%% Building the model
mean = 3 * df.inter_miss + 1.3 * df.pc_1 + 0.7 * df.pc_2 -0.5*df.pc_3 -0.535 * df.tidyness

# Create dummy for site
df.site = np.array(df.site, dtype = str)

sites = np.array(pd.get_dummies(df.site))

#%%
sg = 0.8
ss = 0.1
se = 0.1
S = np.matmul(sites, sites.T)

V = sg * A  + ss * S + se * np.identity(5000)

#%% simulate random effects
y = np.random.multivariate_normal(mean = mean, cov = V)

