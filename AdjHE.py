#! /usr/bin/env python3

"""
AdjHE estimator
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-06-06
"""

##############################################################
# The main file of the AdjHE estimator: loads, cleans, selects
# , loops, and store heritability estimates 
##############################################################

#%% For troubleshooting 
import os
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
c_args= {}
c_args['argfile'] = "Example/Argfile.json"

#%%

import numpy as np
import pandas as pd
import timeit
import itertools
from functions.AdjHE_estimator import load_n_estimate
from functions.load_data import ReadGRMBin, multirange, load_data
from functions.parser import get_args, read_flags
from functions.traits_visualizer import covs_vs_cov_of_interest
# from pathlib import Path


#%%
# Get command line arguments
c_args = get_args()
# convert the arguments to usable Python objects in a dictionary
args = read_flags(c_args)
print(args)

# Save each dictionary item as it's own object to make it easier to reference in the code below
prefix = args["prefix"]
covar = args["covar"]
pheno = args["pheno"]
mpheno = args["mpheno"]
PC = args["PC"]
npc = args["npc"]
out = args["out"]
std = args["std"]
k = args["k"]
covars = args["covars"]
loop_covs = args["loop_covars"]
PredLMM = args["PredLMM"]
RV = args["RV"]


# %% Read GRM
print("Reading GRM: ", prefix)

# Time reading the GRM and other data
start_read = timeit.default_timer()

# Read in grm
G = ReadGRMBin(prefix)
# Get specific detials about the GRM
ids = G['id']
n_phen_nona = G['n_phen_nona']
GRM_array_nona = np.zeros((n_phen_nona, n_phen_nona))
GRM_array_nona[np.diag_indices(n_phen_nona)] = G['diag']

###############################
# Don't know what this is doing
if(k == 0):
    k = n_phen_nona
temp_i = 0
temp = 0
# k= args.k

l = list(range(k, n_phen_nona, k))
l.append(n_phen_nona)
for i in l:
    cor = multirange(range(temp_i, i))
    GRM_array_nona[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
    GRM_array_nona.T[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
    temp = temp + len(cor['b'])
    del(cor)
    temp_i = i
################################


df, covariates, phenotypes = load_data(pheno_file=pheno, cov_file=covar, PC_file=PC)
end_read = timeit.default_timer()
read_time = end_read - start_read

print("It took " + str(read_time) + " (s) to read GRM, covariates, and phenotypes")
print("Phenos + Covars:", df.columns)
#%%
# Save images of covariate relations
covs_vs_cov_of_interest(df, RV, covars, out)

#%%
print("Calculating heritibility")

# create empty list to store heritability estimates
results = pd.DataFrame()

covars = [covar-1 for covar in covars]
#%%
# Create the sets of covarates over which we can loop
cov_combos = [covars[0:idx+1] for idx, c in enumerate(covars)]
cov_combos = [list(covariates[cov_combo]) for cov_combo in cov_combos]

# If we don't want to loop, just grab the last item of the generated list assuming the user wants all of those variables included 
if (loop_covs != True) : 
    cov_combos = [cov_combos[-1]]

# get list of phenotype names to regress
mpheno = [phenotypes[i-1] for i in mpheno]
#%%
# loop over all combinations of pcs and phenotypes
for idx, (mp, nnpc, covs) in enumerate(itertools.product(mpheno, npc, cov_combos)):
    r = load_n_estimate(
        df=df, covars=covs, nnpc=nnpc, mp=mp, ids=ids, GRM_array_nona=GRM_array_nona, std=False)
    results = pd.concat([results, r])
    

# %%
print("Writing results")
results.to_csv(out + "/Results/AdjHE.csv", index=False, na_rep='NA')
