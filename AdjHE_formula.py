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

start_time = timeit.default_timer()
parser = argparse.ArgumentParser(prog='Running adjusted HE regression',description='This program gives estimation in formula fashion.\n Make sure you have enough memory to store GRM matrix in python.',formatter_class=RawTextHelpFormatter)
parser.add_argument('--PC', type=str, help='Read PLINK format covariate file contains the PCs \nPCs should be generated using the same set of individuals in GRM files.\nIf --npc is not specified then all PCs in the file will be used.')
parser.set_defaults(PC="NULL")
parser.add_argument('--npc', type=int, help='Specify the number of PCs to be adjusted')
parser.set_defaults(npc=-9)
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


def sum_n_vec(n):
    out = [int(0)] * n
    for i in range(n):
        out[i] = int(((i + 1) * (i + 2) / 2) - 1)
    return(out)

	
def ReadGRMBin(prefix, AllN = False):
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
    if AllN:
        N = np.fromfile(NFileName, dtype = dt)
    else:
        with open(NFileName, mode='rb') as f:
            record = f.read(entry_size)
            N = unpack(entry_format, record)[0]
            N = int(N)
    i = sum_n_vec(n)
    out = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids,'N':N}
    return(out)

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


def myformula1(A,y,trA=None,trA2=None):
    std_y = (y-np.mean(y))/np.std(y)
    if (trA is None) and (trA2 is None):
        trA = np.sum(np.diag(A))
        trA2 = np.sum(np.multiply(A,A))
    n = A.shape[1]
    yay = np.dot(std_y.T,np.dot(A,std_y))
    yty = np.dot(std_y,std_y)
    if (npc==0):
        denominator = trA2 - 2*trA + n
        nominator = n - trA + yay - yty
    else:
        pc = final_PC
        s = np.diag(np.dot(pc.T,np.dot(A,pc)))
        b = s - 1
        c = np.dot(std_y,pc)**2 - 1
        denominator = trA2 - 2*trA + n - np.sum(b**2)
        nominator = n - trA + yay - yty - np.sum(b*c)
    h2 = nominator/denominator
    var_ge = 2/denominator
#    tau = n/nmarkers
#    b1 = (1-np.sqrt(tau))**2
#    b2 = (1+np.sqrt(tau))**2
#    r = b2-b1
#    a1 = h2-1
#    a2 = 1-2*h2
#    trace_A2_MP = 0.5*(r+2*b1)*n
#    trace_A3_MP = (5/16*r**2+b1*b2)*n
#    trace_A4_MP = (7*r**3+30*b1*r**2+48*b1**2*r+32*b1**3)/32*n
#    if (npc==0):
#    #    var_MP = 2/denominator
#        var_ge = 2/denominator
#    else:
#        trace_A_MP = trA - np.sum(s)
#        a = denominator
#    #    var_MP=2/a**2*(h2**2*trace_A4_MP+(n-npc)*a1**2+(a2**2+2*h2*a1)*trace_A2_MP+2*a1*a2*trace_A_MP+2*h2*a2*trace_A3_MP)
#        var_ge = 2/a
    return h2,np.sqrt(var_ge)


def myformula2(A,y,trA=None,trA2=None):
#    y = y - np.mean(y)
    if (trA is None) and (trA2 is None):
        trA = np.sum(np.diag(A))
        trA2 = np.sum(np.multiply(A,A))
    n = A.shape[1]
    yay = np.dot(y,np.dot(A,y))
    yty = np.dot(y,y)
    tn = np.sum(y)**2/n # all 1s PC
    if (npc==0):
        sigg = n*yay - trA*yty
        sigg = sigg-yay+tn*trA # add 1's
        sige = trA2*yty - trA*yay
        sige = sige-tn*trA2 # add 1's
        denominator = trA2 - 2*trA + n
    else:
        pc = final_PC
        pcA = np.dot(pc.T,A)
        pcApc = np.dot(pcA,pc)
        s = np.diag(pcApc) #pciApci
        b = s-1
        t = np.dot(y,pc)**2 #ypcipciy
        a11 = trA2 - np.sum(s**2) 
        a12 = trA - np.sum(s)
        b1 = yay - np.sum(s*t)
        b2 = yty - np.sum(t)
        sigg = (n-npc)*b1 - a12*b2
        sigg = sigg-yay+tn*a12 # add 1's
        sige = a11*b2 - a12*b1
        sige = sige-tn*a11 # add 1's
        denominator = trA2 - 2*trA + n - np.sum(b**2)
    h2 = sigg/(sigg+sige)
    var_ge = 2/denominator
    return h2,np.sqrt(var_ge)

def regout(y):
    X = cov_selected
    XTX_inv = np.linalg.inv(np.dot(X.T,X))
    XTY = np.dot(X.T,y)
    beta = np.dot(XTX_inv,XTY)
    res = y - np.dot(X,beta)
    return(res)


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
phenotypes = pd.DataFrame(np.loadtxt(args.pheno))
phenotypes.index = phenotypes.iloc[:,0].astype("int32")
intersection_indiv = np.intersect1d(ids.iloc[:,0].astype("int32"), phenotypes.iloc[:,0].astype("int32"))
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



if (args.std == True):
    h2,se = myformula1(GRM_array_nona,res_y,trA=trace_A,trA2=trace_A2)
else:
    h2,se= myformula2(GRM_array_nona,res_y,trA=trace_A,trA2=trace_A2)
logging.info('h2: '+str(h2))
logging.info('Standard error: '+str(se))
logging.info('It takes: '+str(timeit.default_timer() - start_time)+' seconds.')
logging.info('Memory usage:'+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
