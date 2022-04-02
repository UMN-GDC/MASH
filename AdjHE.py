import statsmodels.api as sm
import os
import numpy as np
import pandas as pd
import timeit
import resource
from functions.AdjHE_estimator import AdjHE_estimator
#os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas, load_data
from functions.AdjHE_parser import args 
import itertools

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
# npc=[4,5,6,7,8,9,10]
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
df = load_data(pheno_file = pheno, IDs= ids, cov_file=covar, PC_file=PC, npc = npc, covars = covars)
print(covars)
print("Covariates:", df.columns)
print("Calculating heritibility")

# Empty vectors of heritability SEs, time and memory
h2s =np.empty; SEs = []; Mems = []; Times = []
# create empty list to store heritability estimates
results = pd.DataFrame(np.zeros((len(mpheno) * len(npc), 6)))
results.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage"]

#%%
# loop over all combinations of pcs and phenotypes
for idx, (mp, nnpc)in enumerate(itertools.product(mpheno, npc)):
    # Get indices for ID variables
    id_cols = (df.columns == "FID") + (df.columns == "IID") 
    print(idx, mp, nnpc)    
    # Get the full range of pc columns
    pc_cols = ["PC_" + str(p) for p in range(1, nnpc +1)]
    pc_cols = [c in pc_cols for c in df.columns]
    # grab the covariate columns
    covar_cols = df.columns.str.startswith('Covar')
    pheno_col = df.columns == 'Pheno_'+ str(mp)
    # Combine boolean vectors for all selected columns
    all_columns = id_cols + pc_cols + covar_cols + pheno_col
    # select the dataframe
    temp = df.loc[:,all_columns]
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
    results.iloc[idx,0], results.iloc[idx,1] = AdjHE_estimator(A= GRM_nonmissing, data = temp, mp = mp, npc=nnpc, std=std)
    results.iloc[idx, 2] = mp
    results.iloc[idx, 3] = nnpc
    # Get time for each estimate
    results.iloc[idx, 4] = timeit.default_timer() - start_est
    # Get memory for each step (This is a little sketchy)
    results.iloc[idx, 5] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(results.iloc[idx,0])
# %%
# print("Heritability estimate: " + str(h2[0]))
# print("With Standard error: " + str(se))
print("Writing results")
results.to_csv(out + ".csv", index = False, na_rep='NA')
