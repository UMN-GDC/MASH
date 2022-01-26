# Delete this
import statsmodels.api as sm
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas, load_data
import os
import numpy as np
import pandas as pd
import timeit
import resource
from functions.AdjHE_estimator import AdjHE_estimator
#os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.AdjHE_parser import prefix, npc
from functions.AdjHE_parser import *

# good stuff
# from argparse import RawTextHelpFormatter
print(args)
# %% for troubleshooting
#os.chdir("/home/christian/Scripts/Basu_herit")
#prefix = "Example/grm"
#pheno = "Example/pheno.phen"
#covar = "Example/covar.csv"
#PC = "Example/pcas.eigenvec"
#k = 0
#npc = 2
#mpheno = 1
#std = False
#out = "Example/results"

print("reaading GRM")


# %% Read GRM
G = ReadGRMBin(prefix)
ids = G['id']
ids = ids.rename(columns = {0:"FID", 1:"IID"})
n_phen_nona = ids.shape[0]

print("loading data")
#%%
df = load_data(pheno_file = pheno, cov_file=covar, PC_file=PC, npc = npc)
#%%

# only regress out covariates if they are entered
res_y = sm.OLS(endog=df.loc[:,"Pheno_" + str(mpheno)], exog=df.drop("Pheno_" + str(mpheno), 1)).fit().resid

# %%

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
print("Calculating heritibility")
res_y.name = "Residual"

df["Residual"] = res_y

#%%
h2, se = AdjHE_estimator(A=GRM_array_nona, data = df, npc=npc, std=std)
# %%

results = {"h2" : h2,
      "SE" : se,
      "Time for analysis(s)" : timeit.default_timer() - start_read,
      "Memory usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
print(h2)
print(se)

results= pd.DataFrame(results, index =[0])
results.to_csv(out + ".csv", index = False )
