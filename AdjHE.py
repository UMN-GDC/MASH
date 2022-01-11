# Delete this
import statsmodels.api as sm
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas
import os
import numpy as np
import pandas as pd
import timeit
import resource
os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.arg_parser import prefix, npc
from functions.arg_parser import *

# good stuff
# from argparse import RawTextHelpFormatter

# %% for troubleshooting
os.chdir("/home/christian/Research/Stat_gen/AdjHE")
prefix = "Example/grm"
pheno = "Example/pheno.phen"
covar = "Example/covar.csv"
PC = "Example/pcas.eigenvec"
k = 0
npc = 2
mpheno = 1
std = False
out = "Example/results"

# %% Read GRM
G = ReadGRMBin(prefix)
ids = G['id']
ids = ids.rename(columns = {0:"FID", 1:"IID"})
n_phen_nona = ids.shape[0]

# %%  load phenotypes and covariates
# load phenotypes
y = read_datas(pheno)
# %%
# read in covariates if nonnull
try:
    cov_selected = read_datas(covar)
except:
    print("No covariates file specified or specified file is not found or cannot be loaded.")
#%%

# onlyt load pcs if non null
try:
    PCs = read_datas(PC)
    if (npc == -9):
        npc = PCs.shape[1] - 2
    # prune it to only the number of pc's wanted
    PCs.iloc[:, list(range(npc + 2))]

except:
    print("No PC file specified or specified file is not found or cannot be loaded.")

# join PC's and covariates
cov_selected = pd.merge(cov_selected, PCs, on = ["FID", "IID"])
cov_selected = pd.merge(cov_selected, ids, on = ["FID", "IID"])

# only regress out covariates if they are entered
res_y = sm.OLS(endog=y, exog=cov_selected).fit().resid
# this is reshaping 1-D array to vector in numpy, this might cause problems for multivariate regression
# res_y = np.reshape(np.asarray(res_y), -1)

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

h2, se = AdjHE_estimator(A=GRM_array_nona, y=res_y, npc=1, std=std, PCs = PCs)
# %%

results = {"h2" : h2[2],
      "SE" : se,
      "Time for analysis(s)" : timeit.default_timer() - start_read,
      "Memory usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
results= pd.DataFrame(results, index =[0])
results.to_csv(out + ".csv", index = False )
