import os
import numpy as np
import pandas as pd
import timeit
import itertools
from functions.AdjHE_estimator import load_n_estimate
#os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.load_data import ReadGRMBin, multirange, load_data
from functions.AdjHE_parser import args 

#%%

prefix= args["prefix"]
covar= args["covar"]
pheno= args["pheno"]
mpheno= args["mpheno"]
PC= args["PC"]
npc=args["npc"]
out=args["out"]
std = args["std"]
k=args["k"]
covars = args["covars"]

print(prefix)
print("Reading GRM")

#%%
# prefix="Example/grm"
# covar="Example/covar.csv"
# pheno="Example/pheno.phen"
# mpheno=[1, 2, 3]
# PC="Example/pcas.eigenvec"
# npc=[1,2, 4,5,6, 8,10]
# out="delete"
# std=False
# k=0
# covars=[1]


# %% Read GRM
G = ReadGRMBin(prefix)
ids = G['id']
ids = ids.rename(columns = {0:"FID", 1:"IID"})
ids = ids.dropna()
ids["FID"] = ids.FID.astype(int)
n_phen_nona = ids.shape[0]

start_read = timeit.default_timer()
nmarkers = G['N']
x = G['diag'].astype('float64')
n_phen_nona = G['diag'].size
GRM_array_nona = np.zeros((n_phen_nona, n_phen_nona))
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

GRM_array_nona[np.diag_indices(n_phen_nona)] = G['diag']
#%%



print("loading data")
df = load_data(pheno_file = pheno, IDs= ids, cov_file=covar, PC_file=PC)
print(covars)
print("Covariates:", df.columns)
print("Calculating heritibility")

# create empty list to store heritability estimates
results = pd.DataFrame(np.zeros((len(mpheno) * len(npc), 7)))
results.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]

#%%
# loop over all combinations of pcs and phenotypes
for idx, (mp, nnpc)in enumerate(itertools.product(mpheno, npc)):
    results.iloc[idx,:] = load_n_estimate(df = df, covars =covars, nnpc=nnpc, mp=mp, ids=ids, GRM_array_nona=GRM_array_nona, std = False)

#%%
print("Writing results")
results.to_csv(out + ".csv", index = False, na_rep='NA')
