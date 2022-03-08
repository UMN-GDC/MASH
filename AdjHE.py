import statsmodels.api as sm
import os
import numpy as np
import pandas as pd
import timeit
import resource
from functions.AdjHE_estimator import AdjHE_estimator
#os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas, load_data
from functions.AdjHE_parser import prefix, npc, covars
from functions.AdjHE_parser import *

# from argparse import RawTextHelpFormatter
print(args)

#################################################
# %% for troubleshooting Basic example
# os.chdir("/home/christian/Research/Stat_gen/Basu_herit")
#prefix = "Example/grm"
#pheno = "Example/pheno2.phen"
#covar = "Example/covar.csv"
#PC = "Example/pcas.eigenvec"
#k = 0
#npc = 2
#mpheno = [1,2,3]
#std = False
#out = "Example/results"
#covars = [1,2]
###############################
# TROUBLESHOOT ABCD data
# prefix= "/panfs/roc/groups/3/rando149/coffm049/ABCD/workflow/01_Gene_QC/filters/filter1/GRMs/full"
#covar= "/panfs/roc/groups/3/rando149/coffm049/ABCD/workflow/02_Phenotypes/empty_covar.csv"
#pheno= "/panfs/roc/groups/3/rando149/coffm049/ABCD/workflow/02_Phenotypes/Pconns/conn_files_short.files"
#mpheno= [1, 2, 3]
#PC= "/panfs/roc/groups/3/rando149/coffm049/ABCD/workflow/01_Gene_QC/filters/filter1/Eigens/full.eigenvec"
#npc=4
#out="AdjHE_first"
#std = False
#k=0
#ids= "/panfs/roc/groups/3/rando149/coffm049/ABCD/workflow/IDs.txt"
####################################################################

print("reaading GRM")


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
df = load_data(pheno_file = pheno, IDs= ids, cov_file=covar, PC_file=PC, npc = npc, covars = covars)
print(covars)
print("Covariates:", df.columns)
print("Calculating heritibility")

# Empty vectors of heritability SEs, time and memory
h2s =np.empty; SEs = []; Mems = []; Times = []
# create empty list to store heritability estimates
results = pd.DataFrame(np.zeros((len(mpheno), 4)))
results.columns = ["h2", "SE", "Time for analysis(s)", "Memory Usage"]

#%%

for idx, mp in enumerate(mpheno):
    # Save temp with just the phenotype we need (I'm sure this can be written given the hints that python returns
    temp=df.loc[:,(df.columns == "FID") + (df.columns == "IID") + df.columns.str.startswith('PC')+ df.columns.str.startswith('Covar')+ (df.columns == 'Pheno_'+ str(mp))]
    # drop missing values from both phenos and covariates
    temp = temp.dropna()    
    # Save residuals of selected phenotype after regressing out PCs and covars
    temp["res" + str(mp)] = sm.OLS(endog= temp.loc[:,"Pheno_" + str(mp)], exog= temp.loc[:,temp.columns.str.startswith('PC')+ temp.columns.str.startswith('Covar')]).fit().resid
    # keep portion of GRM without missingess
    nonmissing = ids[ids.IID.isin(temp.IID)].index
    GRM_nonmissing = GRM_array_nona[nonmissing,:][:,nonmissing]
    # resutls from mp pheno
    start_est = timeit.default_timer()
    # Get heritability and SE estimates
    results.iloc[idx,0], results.iloc[idx,1] = AdjHE_estimator(A= GRM_nonmissing, data = temp, mp = mp, npc=npc, std=std)
    # Get time for each estimate
    results.iloc[idx, 2] = timeit.default_timer() - start_est
    # Get memory for each step (This is a little sketchy)
    results.iloc[idx, 3] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(results.iloc[idx,0])
# %%
# print("Heritability estimate: " + str(h2[0]))
# print("With Standard error: " + str(se))
print("Writing results")
results.to_csv(out + ".csv", index = False, na_rep='NA')
