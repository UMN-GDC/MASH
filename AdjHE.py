import os
import numpy as np
import pandas as pd
import timeit
import itertools
from functions.AdjHE_estimator import load_n_estimate
# os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.load_data import ReadGRMBin, multirange, load_data
from functions.AdjHE_parser import args
from pathlib import Path


# %%

# Get arguments from the argparser
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

print(prefix)
print("Reading GRM")

# %%
# prefix="Example/grm"
# covar="Example/covar.csv"
# pheno="Example/pheno.phen"
# mpheno=[1, 2, 3]
# PC="Example/pcas.eigenvec"
# npc=[1,2, 4,5,6, 8,10]
# out="delete"
# std=False
# k=0
# covars=[2, 1]

# %% Read GRM
# Time reading the GRM and other data
start_read = timeit.default_timer()

# Read in grm
G = ReadGRMBin(prefix)
# Get specific detials about the GRM
ids = G['id']
n_phen_nona = G['n_phen_nona']
GRM_array_nona = np.zeros((n_phen_nona, n_phen_nona))
GRM_array_nona[np.diag_indices(n_phen_nona)] = G['diag']

temp_i = 0
temp = 0
# k= args.k


if(k == 0):
    k = n_phen_nona

l = list(range(k, n_phen_nona, k))
l.append(n_phen_nona)
for i in l:
    cor = multirange(range(temp_i, i))
    GRM_array_nona[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
    GRM_array_nona.T[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
    temp = temp + len(cor['b'])
    del(cor)
    temp_i = i



df, covariates, phenotypes = load_data(pheno_file=pheno, IDs=ids, cov_file=covar, PC_file=PC)
end_read = timeit.default_timer()
read_time = end_read - start_read

print("It took " + str(read_time) + " (s) to read GRM, covariates, and phenotypes")
print("Phenos + Covars:", df.columns)
print("Calculating heritibility")

# create empty list to store heritability estimates
results = pd.DataFrame()

# Create the sets of covarates over which we can loop
cov_combos = [covars[0:idx+1] for idx, c in enumerate(covars)]
cov_combos = [list(covariates[cov_combo]) for cov_combo in cov_combos]
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
results.to_csv(out + ".csv", index=False, na_rep='NA')
