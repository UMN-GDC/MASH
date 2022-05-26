#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 10:35:56 2022

@author: christian
"""
import itertools
import numpy as np
import pandas as pd 

#%% Inputs
prefix="Example/grm"
covar="Example/covar.txt"
pheno="Example/pheno.phen"
mpheno=[1, 2]
npc=[5]
PC="Example/pcas.eigenvec"
out="Example/AdjHE_results"
std=False
k=0
covars=[2,1]
loop_covs=False

#%%
# read data for the column names (only read the first line)
covs = pd.read_table(covar)

#%%
# Build arg files to paralellize jobs 

covars =np.array(covar)-1)

# Create the sets of covarates over which we can loop
cov_combos = [covars[0:idx+1] for idx, c in enumerate(covars)]
cov_combos = [list(covariates[cov_combo]) for cov_combo in cov_combos]
#%%

# first assume that we're not looping over PC's and covarites This means there needs to be only one npc
if (loop_covs == False and len(npc) == 1) :
    print("we can do simple loop") 
#%%
iis = list(itertools.product(mpheno, npc, [covars]))


#  Create all combinations 
print(iis)


