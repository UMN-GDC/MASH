import pandas as pd
import numpy as np
import sys
import os
import psutil
import scipy.sparse as sp
import scipy.sparse.linalg
import inspect
from scipy.sparse import csr_matrix, rand
from scipy.misc import imsave
from struct import unpack, calcsize
from numpy.linalg import inv
from numpy.linalg import multi_dot
import timeit
import logging
import resource
import argparse
from argparse import RawTextHelpFormatter
from functions import *


start_time = timeit.default_timer()
parser = argparse.ArgumentParser(prog='Running adjusted HE regression',description='This program gives estimation in formula fashion.\n Make sure you have enough memory to store GRM matrix in python.',formatter_class=RawTextHelpFormatter)
parser.add_argument('--PC', type=str, help='Read PLINK format covariate file contains the PCs \nPCs should be generated using the same set of individuals in GRM files.\nIf --npc is not specified then all PCs in the file will be used.')
parser.set_defaults(PC="NULL")
parser.add_argument('--npc', type=int, help='Specify the number of PCs to be adjusted')
parser.set_defaults(npc=0)
parser.add_argument('--prefix', type=str, help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id [required]',required=True)
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')
parser.set_defaults(covar="NULL")
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
parser.add_argument('--mpheno',type=int, default=1,help='Specify which phenotype to use from phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
parser.add_argument('--k',type=int,help='Specify the number of rows in restoring the GRM each time.\n This could affect the computation time and memory especially when sample size is large. If not provide, it will process the whole GRM at one time.')
parser.set_defaults(k=0)
parser.add_argument('--out',type=str, help='Specify the output file name. [required]',required=True)
parser.add_argument('--std',action='store_true',default=False,help='Run SAdj-HE (i.e., with standardization)')

args = parser.parse_args()


if (args.npc == "NULL"):
    npc=0
else:
    npc = args.npc
outprefix = args.out
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG,filename=outprefix+'.log',filemode='a')
for arg, value in sorted(vars(args).items()):
    logging.info("Argument %s: %r", arg, value)

prefix = args.prefix

G = ReadGRMBin(prefix)
ids = G['id']
n_phen_nona = ids.shape[0]
phenotypes = pd.read_csv(args.pheno, sep=" ", header =None  )
# remove null only to make it work because I have a null column, need to figure out how to fix this permantetly later
# phenotypes = phenotypes.iloc[: , :-1]

phenotypes.index = phenotypes.iloc[:,0].astype("int32")
intersection_indiv = np.intersect1d(ids.iloc[:,0].astype("int32"), phenotypes.iloc[:,0].astype("int32"))
final_phen = phenotypes.loc[intersection_indiv]

cov_selected = np.ones(n_phen_nona)

# only load covaraites if nonnull 
if (args.covar!="NULL"):
    covariates = pd.read_csv(args.covar, sep=" ", header = 0)
    covariates.index = covariates.iloc[:,0].astype("int32")
    final_covar = covariates.loc[intersection_indiv]
    final_covar = final_covar.values[:,2:]
    cov_selected = np.column_stack((cov_selected,final_covar))

# onlyt load pcs if non null
if (args.PC != "NULL"):
    PCs = pd.DataFrame(np.loadtxt(args.PC))
    PCs.index = PCs.iloc[:,0].astype("int32")
    if (args.npc == -9):
        npc = PCs.shape[1] - 2
    if (npc != 0):
        final_PC = PCs.loc[intersection_indiv]
        final_PC = final_PC.values[:,2:(2+npc)]
        cov_selected = np.column_stack((cov_selected,final_PC))

y = final_phen.iloc[:,args.mpheno+1]


# only regress out covariates if they are entered
res_y = regout(cov_selected, y, args.covar, args.PC)

start_read = timeit.default_timer()
G = ReadGRMBin(prefix)
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

if (args.std == True):
    h2,se = myformula1(A=GRM_array_nona,y=res_y,trA=trace_A,trA2=trace_A2, npc = npc)
else:
    h2,se= myformula2(A=GRM_array_nona,y=res_y,trA=trace_A,trA2=trace_A2, npc =npc)
logging.info('h2: '+str(h2))
logging.info('Standard error: '+str(se))
logging.info('It takes: '+str(timeit.default_timer() - start_time)+' seconds.')
logging.info('Memory usage:'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
