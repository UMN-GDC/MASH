import pandas as pd
import numpy as np
import sys
import os
import psutil
import scipy.sparse as sp
import scipy.sparse.linalg
import math as math
from scipy.misc import imsave
from struct import unpack, calcsize
from numpy.linalg import inv
import logging
import timeit
import argparse
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(prog='Running adjusted HE regression for large samplesize in parallel (2)',description='IMPORTANT: Run this for heritability estimation only after finishing HE_reg_s1.py',formatter_class=RawTextHelpFormatter)
parser.add_argument('--Npart',type=int, help='Specify the total number of part-GRM files. [required]\n e.g. 200',required=True)
parser.add_argument('--id', type=str, help='Read GCTA format .grm.id file contains all individuals. You can generate this file by "cat partgrm_prefix_*.grm.id > all_individuals.grm.id" [required]',required=True)
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
parser.add_argument('--mpheno',type=int, default=1,help='Specify which phenotype to use from phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
parser.add_argument('--out',type=str, help='Specify the output directory name. This should be the SAME as the one used in HE_reg_s1.py',required=True)

parser.add_argument('--PC', type=str, help='Read PLINK format covariate file contains the PCs \nPCs should be generated using the same set of individuals in GRM files. ID should be matched with --id \nIf --npc is not specified then all PCs in the file will be used.')
parser.set_defaults(PC="NULL")
parser.add_argument('--npc', type=int, help='Specify the number of PCs to be adjusted')
parser.set_defaults(npc=-9)
parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates BESIDES PCs to be adjusted')
parser.set_defaults(covar="NULL")
parser.add_argument('--std',action='store_true',default=False,help='Run SAdj-HE regression (i.e., with standardization)')



args = parser.parse_args()
start_time0 = timeit.default_timer()

def multirange(counts):
    counts = np.asarray(counts)
    # Remove the following line if counts is always strictly positive.
    counts = counts[counts != 0]
    counts1 = counts[:-1]
    reset_index = np.cumsum(counts1)
    incr = np.ones(counts.sum(), dtype=int)
    incr[0] = 0 
    incr[reset_index] = 1 - counts1
    # Reuse the incr array for the final result.
    incr.cumsum(out=incr)
    out = {'a':incr,'b':np.repeat(counts,counts)}
    return(out)

def sum_n_vec(m,n):
    out = [int(0)] * n
    for i in range(n):
        out[i] = int(((i + 1) * (i + 2) / 2) - 1+ (m) * (i + 1))
    return(out)
    
        
def ReadGRMBin(prefix, m=0,AllN = False):
    BinFileName  = prefix + ".grm.bin"
    NFileName = prefix + ".grm.N.bin"
    IDFileName = prefix + ".grm.id"
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
    entry_format = 'f' # N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Read IDs
    ids = pd.read_csv(IDFileName, sep = '\t', header = None)
    ids_vec = ids.iloc[:,1]
    n = len(ids.index)
    ids_diag = ['NA' for x in range(n)]
    n_off = int(n * (n - 1) / 2) 
    ## Read relatedness values
    grm = np.fromfile(BinFileName, dtype = dt)
    ## Read number of markers values
#    if AllN:
#        N = np.fromfile(NFileName, dtype = dt)
#    else:
#        with open(NFileName, mode='rb') as f:
#            record = f.read(entry_size)
#            N = unpack(entry_format, record)[0]
#            N = int(N)
    i = sum_n_vec(m,n)
    out = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids}
    return(out)
    
def smartway(x1,x2):
    x12 = x1*x2
    temp_sum = np.sum(x12)
    lower_sum = 0
    for i in range(len(x12)):
        temp_sum = temp_sum - x12[i]
        lower_sum = lower_sum + x12[i] * temp_sum
    return(lower_sum)

def outindex(counter,temp_n):
    temp_count = int(temp_n*(temp_n-1)/2)
    temp_array1 = np.array(range(temp_count))
    temp_array2 = np.array([counter*(i+2) for i in range(temp_n-1)])
    temp_array2 = np.repeat(temp_array2,np.array(range(temp_n-1))+1)
    out_index = temp_array1+temp_array2
    return(out_index)

def fun1(current_sum,counter,n,pc,lower_diag):
    temp_out = outindex(counter,n)
    temp_pc = pc[counter:counter+n]
    pre_pc = pc[:counter]
    temp_outer1 = np.outer(temp_pc,temp_pc)
    temp_outer2 = np.outer(temp_pc,pre_pc)
    cor = multirange(range(n))
    temp_sum1 = np.dot(lower_diag[temp_out],temp_outer1[cor['b'],cor['a']])
    temp_sum2 = np.dot(np.delete(lower_diag,temp_out),temp_outer2.ravel())
    temp_sum = temp_sum1 + temp_sum2
    current_sum = current_sum + temp_sum
    return(current_sum)
        

def regout(y):
    X = cov_selected
    XTX_inv = np.linalg.inv(np.dot(X.T,X))
    XTY = np.dot(X.T,y)
    beta = np.dot(XTX_inv,XTY)
    res = y - np.dot(X,beta)
    return(res)

outprefix = args.out
npc = args.npc
os.chdir(outprefix)
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG,filename='pheno'+str(args.mpheno)+'.log',filemode='a')
#logging.info('begin..!')

ids = pd.read_csv(args.id,sep='\s+',header=None)
n_phen_nona = ids.shape[0]
phenotypes = pd.DataFrame(pd.read_csv(args.pheno,sep='\s+',header=None))
final_phen = pd.merge(ids,phenotypes,how='inner',on=[0,1])


#logging.info('file read.. finish!')
cov_selected = np.ones(n_phen_nona)
if (args.covar!="NULL"):
    covariates = pd.DataFrame(pd.read_csv(args.covar,sep='\s+',header=None))
    final_covar = pd.merge(final_phen,covariates,how='inner',on=[0,1])
    final_covar = final_covar.values[:,2:]
    cov_selected = np.column_stack((cov_selected,final_covar))
if (args.PC != "NULL"):
    PCs = pd.DataFrame(pd.read_csv(args.PC,sep='\s+',header=None))
    final_PC = pd.merge(final_phen,PCs,how='inner',on=[0,1])
    if (args.npc == -9):
        npc = PCs.shape[1] - 2
    if (npc != 0):
        final_PC = PCs.loc[intersection_indiv]
        final_PC = final_PC.values[:,2:(2+npc)]
        cov_selected = np.column_stack((cov_selected,final_PC))

fun1_sum = 0

for igrm in range(1,args.Npart+1):
    f = 'pheno'+str(args.mpheno)+'_'+str(igrm)
    curr = np.loadtxt(f)
    fun1_sum = fun1_sum + curr

y = final_phen.values[:,args.mpheno+1]
res_y = regout(y)
std_y = (res_y - np.mean(res_y))/np.std(res_y)

if(args.std==False):
    xtx = np.zeros((npc+2,npc+2))
    xty = np.zeros(npc+2)
    ##calculate covariate&covariate
    for i in range(npc):
        temp_xty = smartway(final_PC[:,i],res_y)
        xty[i] = temp_xty + np.dot(final_PC[:,i]**2,res_y**2) 
        xtx[i,-2] = np.sum(final_PC[:,i]**2)
        for j in range(i,npc):
            temp_xtx = smartway(final_PC[:,i],final_PC[:,j])
            xtx[i,j] = temp_xtx + np.dot(final_PC[:,i]**2,final_PC[:,j]**2)
    
    xty[-2] = np.sum(res_y**2)
    xtx[-2,-2] = y.shape[0]
    
    logging.info('PCs-PCs compute done; PCs-y comupute done. It takes: '+str(timeit.default_timer() - start_time0)+' seconds.')
    
    #np.savetxt(outprefix+'.xtx',xtx)
    #np.savetxt(outprefix+'.xty',xty)
    
    xty[-1] = fun1_sum[-3] #grm_y
    xtx[-1,-1] = fun1_sum[-1] #grm_g
    xtx[-2,-1] = fun1_sum[-2] #grm_e
    xtx[:-2,-1] = fun1_sum[:-3] #grm_pc
    
    idx = np.tril_indices(npc+2,-1,npc+2)
    xtx[idx] = xtx.T[idx]
    betas = np.dot(np.linalg.inv(xtx), xty)
    logging.info('h2: '+str(betas[-1]/(betas[-1]+betas[-2])))
    
    o = np.ones(n_phen_nona)
    nstar = (n_phen_nona*(n_phen_nona-1)/2)
    sigma_1 = np.sqrt((smartway(res_y,res_y)-smartway(res_y,o)**2/nstar)/nstar)
    se_1 = np.sqrt(sigma_1 * np.linalg.inv(xtx)[-1,-1])
    sigma_p = np.std(res_y)
    logging.info('standard error: '+str(se_1/sigma_p))
    logging.info('Vg: '+str(betas[-1]))
    logging.info('Ve: '+str(betas[-2]))
    #logging.info('finally done... it takes: '+str(timeit.default_timer() - start_time0)+' seconds.'+'\n')
else:
    xtx = np.zeros((npc+1,npc+1))
    xty = np.zeros(npc+1)
    ##calculate covariate&covariate
    for i in range(npc):
        xty[i] = smartway(final_PC[:,i],std_y)
        for j in range(i,npc):
            xtx[i,j] = smartway(final_PC[:,i],final_PC[:,j])
    
    logging.info('PCs-PCs compute done; PCs-y comupute done. It takes: '+str(timeit.default_timer() - start_time0)+' seconds.')
    
    
    xty[-1] = fun1_sum[-2]
    xtx[-1,-1] = fun1_sum[-1]
    xtx[:-1,-1] = fun1_sum[:-2]
    
    idx = np.tril_indices(npc+1,-1,npc+1)
    xtx[idx] = xtx.T[idx]
    betas = np.dot(np.linalg.inv(xtx), xty)
    logging.info('h2: '+str(betas[-1]))
    
    o = np.ones(n_phen_nona)
    nstar = (n_phen_nona*(n_phen_nona-1)/2)
    sigma_1 = np.sqrt((smartway(std_y,std_y)-smartway(std_y,o)**2/nstar)/nstar)
    se_1 = np.sqrt(sigma_1 * np.linalg.inv(xtx)[-1,-1])
    logging.info('standard error: '+str(se_1))
