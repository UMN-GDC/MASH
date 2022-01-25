#-----------------------Loading  the required Module------------------------------------------
import os
import timeit
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas
from functions.PredLMM_estimator import derivative_minim_sub, derivative_minim_full
import logging
import resource
import argparse
from argparse import RawTextHelpFormatter
from scipy.sparse import csr_matrix, rand
import numpy as np
from numpy.linalg import inv
from scipy.optimize import newton
from scipy.linalg.blas import dgemm,sgemm,sgemv
from copy import copy
import pandas as pd 
from functions.PredLMM_parser import *
from functions.PredLMM_parser import *


# logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
#                     level=logging.DEBUG,filename=outprefix+'.log',filemode='a')
# for arg, value in sorted(vars(args).items()):
#     logging.info("Argument %s: %r", arg, value)

start_read = timeit.default_timer()
#%%
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




y = read_datas(args.pheno) 
X = read_datas(args.covar)

#%%

#----------------------Knot selection and selecting corresponding vectors----------------------------
subsample_size = 500;
sub_sample = sorted(np.random.choice(range(0,N),subsample_size,replace=False))
non_subsample = np.setdiff1d(range(0,N),sub_sample)
indices = np.hstack((sub_sample,non_subsample))
GRM_array = np.float32(GRM_array[np.ix_(indices,indices)].T)
y = y.iloc[indices]; X=X.iloc[indices]; X_T = X.T;
#%%
G_selected = GRM_array[range(0,subsample_size),:][:,range(0,subsample_size)]

phen_sub = np.array(y.iloc[range(0,subsample_size)].drop(["FID", "IID"], axis = 1))
cov_sub = np.array(X.iloc[range(0,subsample_size)].drop(["FID","IID"], axis =1))
#%%
# logging.info('knot selection and corresponding vectors seleected: '+str(timeit.default_timer() - start_read)+' seconds.')
Knot_sel_time = timeit.default_timer() - start_read
#%%
#------------------Fitting LMM using only the selected subsample (set of knots)-------------------------
A_selc = np.copy(G_selected)-np.identity(subsample_size)
result_subsample = derivative_minim_sub(phen_sub, cov_sub, cov_sub.T, G_selected, A_selc, subsample_size)
# print(result_subsample)

#%%
#------------------Running PredLMM----------------------------------------------------------------------
Ct =  np.copy(GRM_array[range(0,subsample_size),:],order='F')
C12 = Ct[:,range(subsample_size,N)]
id_diag = np.diag_indices(N)
diag_G_sub = GRM_array[id_diag]
G_inv = inv(G_selected).T
GRM_array[np.ix_(range(subsample_size,N),range(subsample_size,N))] = sgemm(alpha=1,a=C12.T,b=sgemm(alpha=1,a=G_inv,b=C12))
#%%
del G_inv, C12
#%%
add = copy(-GRM_array[id_diag] + diag_G_sub) ## diagonal adjustment
np.fill_diagonal(GRM_array, - 1 + diag_G_sub)
#%%
Xnew = X.drop(["FID", "IID"], axis =1)
ynew = y.drop(["FID", "IID"], axis =1)
result_full = derivative_minim_full(ynew, Xnew, Xnew.T, Ct, id_diag, add, G_selected, GRM_array, N)
# print(result_full)
#%%
GREML_sub_est = result_subsample['Heritability estimate'][0,0]
GREML_sub_sd = result_subsample['SD of heritability estimate']
GREML_sub_var = result_subsample['Variance estimate'][0][0]
#%%
Pred_est = result_full['Heritability estimate'][0][0]
Pred_sd = result_full['SD of heritability estimate']
Pred_var = result_full['Variance estimate'][0][0]




#%%
sub_results = {"Est" : GREML_sub_est,
               "SD" : GREML_sub_sd,
               "Var" : GREML_sub_var,
               "Knot time" : Knot_sel_time,
               "Time for analysis(s)" : timeit.default_timer() - start_read,
               "Memory usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
results = pd.DataFrame(sub_results, index = ["GREML"])
#%%
full_results = {"Est" : Pred_est, 
                "SD" : Pred_sd,
                "Var" : Pred_var,
                "Knot time" : Knot_sel_time,
                "Time for analysis(s)" : timeit.default_timer() - start_read,
                "Memory usage" : resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
temp_results = pd.DataFrame(full_results, index = ["PredLMM"])
#%%
temp = [results, temp_results]
results = pd.concat(temp)
# write indices because that has estimators name
results.to_csv(args.out + ".csv", index = True )

