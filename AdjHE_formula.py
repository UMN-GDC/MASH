# Delete this 
import os
import numpy as np
import pandas as pd
import timeit
import resource
os.chdir("/home/christian/Research/Stat_gen/AdjHE/")
prefix = "Example/grm" 

# good stuff
import argparse
# from argparse import RawTextHelpFormatter
from functions.load_data import sum_n_vec, ReadGRMBin, multirange, read_datas
from functions.estimate import AdjHE_estimator
from functions.arg_parser import args, prefix, npc, logging
import statsmodels.api as sm


#%% Read GRM
G = ReadGRMBin(prefix)
ids = G['id']
n_phen_nona = ids.shape[0]

#%%
# seed the covariates matrix with a column of 1's for the intercept
cov_selected = np.ones(n_phen_nona)

# load phenotypes and covariates
y = read_datas(args.pheno, ids)

# read in covariates if nonnull
if (args.covar != "NULL"):
 X = read_datas(args.covar, ids)
 # stack the covariates onto the incercepts
 cov_selected = np.column_stack((cov_selected,X))
 
 #%%

# onlyt load pcs if non null
if (args.PC != "NULL"):
    PCs = pd.DataFrame(np.loadtxt(args.PC))
    PCs.index = PCs.iloc[:,0].astype("int32")
    if (args.npc == -9):
        npc = PCs.shape[1] - 2
    if (args.npc != 0):
        final_PC = PCs.loc[intersection_indiv]
        final_PC = final_PC.values[:,2:(2+npc)]
        cov_selected = np.column_stack((cov_selected,final_PC))

# y = y.iloc[:,args.mpheno+1]
#%%
cov_selected = pd.DataFrame(cov_selected)
# only regress out covariates if they are entered
res_y = sm.OLS(endog=y, exog = X).fit().resid
# this is reshaping 1-D array to vector in numpy, this might cause problems for multivariate regression
# res_y = np.reshape(np.asarray(res_y), -1)

#%%

start_read = timeit.default_timer()
nmarkers = G['N']
x = G['diag'].astype('float64')
n_phen_nona = G['diag'].size
GRM_array_nona = np.zeros((n_phen_nona,n_phen_nona))
temp_i = 0
temp  = 0
k = args.k

if(k == 0):
    k = n_phen_nona

l = list(range(k,n_phen_nona,k))
l.append(n_phen_nona)
for i in l:
    cor = multirange(range(temp_i,i))
    GRM_array_nona[cor['b'],cor['a']] = G['off'][temp:temp+len(cor['b'])]
    GRM_array_nona.T[cor['b'],cor['a']] = G['off'][temp:temp+len(cor['b'])]
    temp = temp + len(cor['b'])
    del(cor)
    temp_i = i

GRM_array_nona[np.diag_indices(n_phen_nona)] = G['diag']
logging.info('GRM matrix restored done. It takes: '+str(timeit.default_timer() - start_read)+' seconds.')


trace_A = np.sum(x)
trace_A2 = 2*np.sum(G['off'].astype('float64')**2) + np.sum(x**2)
#print(trace_A)
#print(trace_A2)
del(G)

# print(res_y.head())

h2,se = AdjHE_estimator(A=GRM_array_nona,y=res_y,trA=trace_A,trA2=trace_A2, npc = npc, std = args.std)
logging.info('h2: '+str(h2))
logging.info('Standard error: '+str(se))
logging.info('It takes: '+str(timeit.default_timer() - start_read)+' seconds.')
logging.info('Memory usage:'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
