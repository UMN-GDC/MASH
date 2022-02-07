# Delete this
import statsmodels.api as sm
import os
import numpy as np
import pandas as pd
import timeit
import resource
from functions.AdjHE_estimator import AdjHE_estimator
#os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas, load_data
from functions.AdjHE_parser import prefix, npc
from functions.AdjHE_parser import *

# good stuff
# from argparse import RawTextHelpFormatter
print(args)
# %% for troubleshooting
# os.chdir("/home/christian/Scripts/Basu_herit")
# prefix = "Example/grm"
# pheno = "Example/pheno2.phen"
# covar = "Example/covar.csv"
# PC = "Example/pcas.eigenvec"
# k = 0
# npc = 2
# mpheno = [1,2,3]
# std = False
# out = "Example/results"

print("reaading GRM")


# %% Read GRM
G = ReadGRMBin(prefix)
ids = G['id']
ids = ids.rename(columns = {0:"FID", 1:"IID"})
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
df = load_data(pheno_file = pheno, cov_file=covar, PC_file=PC, npc = npc)
# dropping nas for ease of use
df = df.dropna()

print("Calculating heritibility")

# keep portion of GRM without missingess
nonmissing = ids.IID.isin(df.IID)
GRM_nonmissing = GRM_array_nona[nonmissing,:][:,nonmissing]
#%%
# regress covariates 
for mp in mpheno :
    df["res"+str(mp)] = sm.OLS(endog=df.loc[:,"Pheno_" + str(mp)], 
                               exog=df.loc[df.columns.str.startswith('PC')+ df.columns.str.startswith('Covar')]).fit().resid
# remove nonresidualized phenos
df = df.loc[:, ~df.columns.str.startswith('Pheno')]
#%%
# select everything else besides y
X = df.loc[:, ~df.columns.str.startswith('res')]
# Empty vectors of heritability SEs, time and memory
h2s =np.empty; SEs = []; Mems = []; Times = []
#%%

results = pd.DataFrame(np.zeros((len(mpheno), 4)))
results.columns = ["h2", "SE", "Time for analysis(s)", "Memory Usage"]

#%%

for mp in mpheno :
    # merge temporary dataframe
    tempy = df[["FID", "IID", "res" + str(mp)]]
    # rename to make easier
    tempy.columns = ["FID", "IID", "Residual"]
    temp = pd.merge(X, tempy, on = ["FID", "IID"])
    results.iloc[(mp-1),0], results.iloc[(mp-1),1] = AdjHE_estimator(A= GRM_nonmissing, data = temp, npc=npc, std=std)
    results.iloc[(mp-1), 2] = timeit.default_timer() - start_read
    results.iloc[(mp-1), 3] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(results.iloc[(mp-1),0])
# %%
# print("Heritability estimate: " + str(h2[0]))
# print("With Standard error: " + str(se))
print("Writing results")
results.to_csv(out + ".csv", index = False )
