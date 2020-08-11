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

parser = argparse.ArgumentParser(prog='Running adjusted HE regression for large sample size in parallel (1)',description='IMPORTANT: Run this for heritability estimation before HE_reg_s2.py',epilog='After processing all GRMs then run HE_reg_s2.py to get the estimation',formatter_class=RawTextHelpFormatter)
parser.add_argument('--prefix', type=str, help='prefix for GCTA format part GRM files [required]\n e.g. /path/to/partgrm/ukbiobank.part_200_',required=True)
parser.add_argument('--Npart',type=int, help='Specify the total number of part-GRM files. [required]\n e.g. 200',required=True) 
parser.add_argument('--job', type=int, help='Specify which part of GRM to process this time. Make jobs run in parallel in your system. [required]\n e.g. ${PBS_ARRAYID}',required=True)
parser.add_argument('--id', type=str, help='Read GCTA format .grm.id file contains all individuals [required]',required=True)
parser.add_argument('--pheno', type=str, help='Read PLINK format phenotype file [required]\nIf --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.', required=True)
parser.add_argument('--mpheno',type=int, default=1,help='Specify which phenotype to use from phenotype file (1 phenotype only)')
parser.set_defaults(mpheno=1)
parser.add_argument('--out',type=str, help='Specify the output directory name. It will save the intermediate results for each part-GRM. Please make sure the directory already exists! [required]',required=True)

parser.add_argument('--covar', type=str, help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')
parser.set_defaults(covar="NULL")
parser.add_argument('--PC', type=str, help='Read PLINK format covariate file contains the PCs \nPCs should be generated using the same set of individuals in GRM files. ID should be matched with --id \nIf --npc is not specified then all PCs in the file will be used.')
parser.set_defaults(PC="NULL")
parser.add_argument('--npc', type=int, help='Specify the number of PCs to be adjusted')
parser.set_defaults(npc=-9)
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
npc = args.npc
outprefix = args.out

#if not os.path.exists(outprefix):
#    os.mkdir(outprefix)

os.chdir(outprefix)

logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG,filename='pheno'+str(args.mpheno)+'.log',filemode='a')

ids = np.loadtxt(args.id)
n_phen_nona = ids.shape[0]
phenotypes = pd.DataFrame(np.loadtxt(args.pheno))
phenotypes.index = phenotypes.iloc[:,0].astype("int32")
intersection_indiv = np.intersect1d(ids[:,0].astype("int32"), phenotypes.iloc[:,0].astype("int32"))
final_phen = phenotypes.loc[intersection_indiv]

cov_selected = np.ones(n_phen_nona)
if (args.covar!="NULL"):
    covariates = pd.DataFrame(np.loadtxt(args.covar))
    covariates.index = covariates.iloc[:,0].astype("int32")
    final_covar = covariates.loc[intersection_indiv]
    final_covar = final_covar.values[:,2:]
    cov_selected = np.column_stack((cov_selected,final_covar))
if (args.PC != "NULL"):
    PCs = pd.DataFrame(np.loadtxt(args.PC))
    PCs.index = PCs.iloc[:,0].astype("int32")
    if (args.npc == -9):
        npc = PCs.shape[1] - 2
    if (npc != 0):
        final_PC = PCs.loc[intersection_indiv]
        final_PC = final_PC.values[:,2:(2+npc)]
        cov_selected = np.column_stack((cov_selected,final_PC))


y = final_phen.values[:,args.mpheno+1]
res_y = regout(y)
std_y = (res_y - np.mean(res_y))/np.std(res_y)


current_sum_y = 0
current_sum_grm = 0
current_sum_grm_e = 0
current_sum_pc = np.zeros(npc)
#prefix='/home/saonli/shared/bb/dosage_data3/grmblock/316k_unrelated_566k.part_200_'
prefix = args.prefix


igrm = args.job
iformat = len(str(args.Npart))
exec('prefix1 = prefix + \'{:0'+str(iformat)+'}\'.format(igrm)')
tempid = np.loadtxt(prefix1+'.grm.id')[:,0].astype("int32")
tempid_1 = tempid[0]
counter = int(np.where(ids[:,0] == tempid_1)[0])
G = ReadGRMBin(prefix1, m = counter)
n_ind = len(G['id'])
lower_diag = G['off'].astype('float64')
if (args.wostd==False):
    diag = G['diag'].astype('float64')
    current_sum_grm = np.dot(lower_diag,lower_diag) + np.dot(diag,diag) 
    current_sum_grm_e = np.sum(diag)
    current_sum_y = fun1(current_sum_y,counter,n_ind,res_y,lower_diag) + np.dot(diag,res_y[counter:counter+n_ind]**2)
    for j in range(npc):
        current_sum_pc[j] = fun1(current_sum_pc[j],counter,n_ind,final_PC[:,j],lower_diag) + np.dot(diag,final_PC[counter:counter+n_ind,j]**2)
    current_out = np.append(current_sum_pc,[current_sum_y,current_sum_grm_e,current_sum_grm])
else:
    current_sum_grm = np.dot(lower_diag,lower_diag)
    current_sum_y = fun1(current_sum_y,counter,n_ind,std_y,lower_diag)
    for j in range(npc):
        current_sum_pc[j] = fun1(current_sum_pc[j],counter,n_ind,final_PC[:,j],lower_diag)
    current_out = np.append(current_sum_pc,[current_sum_y,current_sum_grm])
logging.info(str(igrm)+'-th GRM loaded and computed done... It takes: ' + str(timeit.default_timer() - start_time0)+' seconds.')

f = open('pheno'+str(args.mpheno)+'_'+str(igrm),'w')
np.savetxt(f,current_out,newline=" ")
f.write('\n')
f.close()


