#-----------------------Loading  the required Module------------------------------------------
import os
import timeit
import resource
from scipy.sparse import csr_matrix
import numpy as np
from numpy.linalg import inv
from scipy.linalg.blas import sgemm
from copy import copy
import pandas as pd 
from functions.load_data import ReadGRMBin, multirange, load_data
from functions.PredLMM_estimator import derivative_minim_sub, derivative_minim_full
from functions.parser import prefix, covar, pheno, out, PC, args, mpheno



#%%

print("Reading GRM")
start_read = timeit.default_timer()
G = ReadGRMBin(prefix)
N = len(G['diag'])
GRM = csr_matrix((N, N));GRM_array = GRM.todense().A1.reshape(N, N)
idx = np.tril_indices(N,-1,N);idy = np.triu_indices(N,1,N);id_diag = np.diag_indices(N)
GRM_array[idx] = G['off'];GRM_array[id_diag] = G['diag'];GRM_array[idy] = GRM_array.T[idy]
GRM_array = np.float32(GRM_array)

#%%

#-----------------------convert the GRM to h5py format for faster loading------------------- 
#hf = h5py.File('Data/example_grm.h5', 'w')
#hf.create_dataset('dataset_1', data=GRM_array)
#hf.close()

#-----------------------loading GRM in h5py format------------------------------------------- 
#hf = h5py.File('Data/example_grm.h5', 'r')
#GRM_array= np.array(hf.get('GRM'),dtype="float32")


print("Reading covariate and phenotype data")
df, covariates, phenotypes = load_data(pheno_file=pheno, cov_file=covar, PC_file=PC)
#%%
#----------------------Knot selection and selecting corresponding vectors----------------------------
print("selecting knots")
subsample_size = 500;
sub_sample = np.random.choice(range(0,N),subsample_size,replace=False)
df_sub = df.iloc[sub_sample, ]
G_selected = GRM_array[sub_sample,:][:,sub_sample]
Knot_sel_time = timeit.default_timer() - start_read
#%%
#------------------Fitting LMM using only the selected subsample (set of knots)-------------------------
print("fitting subsample")
A_selc = np.copy(G_selected)-np.identity(subsample_size)

#%% 
# Temporary
y_sub = np.array(df_sub[phenotypes[0]])
X_sub = np.array(df_sub[covariates[0]])
#%%

result_subsample = derivative_minim_sub(df_sub, G_selected, A_selc, ["inter"], "head_size")
print(result_subsample)

#%%
#------------------Running PredLMM----------------------------------------------------------------------
print("fitting full")
Ct =  np.copy(GRM_array[range(0,subsample_size),:],order='F')
C12 = Ct[:,range(subsample_size,N)]
id_diag = np.diag_indices(N)
diag_G_sub = GRM_array[id_diag]
G_inv = inv(G_selected).T
GRM_array[np.ix_(range(subsample_size,N),range(subsample_size,N))] = sgemm(alpha=1,a=C12.T,b=sgemm(alpha=1,a=G_inv,b=C12))
#%%
# del G_inv, C12
#%%
add = copy(-GRM_array[id_diag] + diag_G_sub) ## diagonal adjustment
np.fill_diagonal(GRM_array, - 1 + diag_G_sub)
#%%
Xnew = X.drop(["fid", "iid"], axis =1)
ynew = y.drop(["fid", "iid"], axis =1)
#%%
Xnew = df["liver_purity"]
ynew = df["inter_miss"]
#%%

result_full = derivative_minim_full(ynew, Xnew, Xnew.T, Ct, id_diag, add, G_selected, GRM_array, N)

#%%

results = pd.DataFrame(np.zeros((len(mpheno), 4)))
results.columns = ["h2", "SE", "Var", "Time"]

Xnew = X.iloc[:,0]
for mp in mpheno:
    result_full = derivative_minim_full(ynew[phenotypes[mp]], Xnew, Xnew.T, Ct, id_diag, add, G_selected, GRM_array, N)
    results.iloc[(mp-1),0] = result_full["Heritability estimate"][0,0]
    results.iloc[(mp-1),1] = result_full["SD of heritability estimate"]
    results.iloc[(mp-1),2] = result_full["Variance estimate"][0,0]
    results.iloc[(mp-1),3] = result_full["Time taken"]
    print(result_full)

results.to_csv(out + ".csv", index= False)

# #%%
# GREML_sub_est = result_subsample['Heritability estimate'][0,0]
# GREML_sub_sd = result_subsample['SD of heritability estimate']
# GREML_sub_var = result_subsample['Variance estimate'][0][0]
# #%%
# Pred_est = result_full['Heritability estimate'][0][0]
# Pred_sd = result_full['SD of heritability estimate']
# Pred_var = result_full['Variance estimate'][0][0]


# print("writing results")

# #%%
# sub_results = {"Est" : GREML_sub_est,
#                "SD" : GREML_sub_sd,
#                "Var" : GREML_sub_var,
#                "Knot time" : Knot_sel_time,
#                "Time for analysis(s)" : timeit.default_timer() - start_read,
#                "Memory usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
# results = pd.DataFrame(sub_results, index = ["GREML"])
# #%%
# full_results = {"Est" : Pred_est, 
#                 "SD" : Pred_sd,
#                 "Var" : Pred_var,
#                 "Knot time" : Knot_sel_time,
#                 "Time for analysis(s)" : timeit.default_timer() - start_read,
#                 "Memory usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
# temp_results = pd.DataFrame(full_results, index = ["PredLMM"])
# #%%
# temp = [results, temp_results]
# results = pd.concat(temp)
# # write indices because that has estimators name
# results.to_csv(args.out + ".csv", index = True )

