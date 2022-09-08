import pandas as pd
import numpy as np
import sys
import os
import psutil
import scipy.sparse as sp
import scipy.sparse.linalg
import inspect
#from cvxpy import*
from scipy.sparse import csr_matrix, rand
from scipy.misc import imsave
from struct import unpack, calcsize
from numpy.linalg import inv
from numpy.linalg import multi_dot
from scipy.linalg import inv
from random import sample 

import timeit
start_time = timeit.default_timer()
os.chdir(path)

covar = sys.argv[1]
prefix = str(sys.argv[3])
pheno = sys.argv[2]
output = str(sys.argv[4])
n_covar = int(sys.argv[5])

phenotypes = pd.DataFrame(np.loadtxt(str(pheno)))
covariates = pd.DataFrame(np.loadtxt(str(covar)))
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
   # if AllN:
   #     N = np.fromfile(NFileName, dtype = dt)
   # else:
   #     with open(NFileName, mode='rb') as f:
   #         record = f.read(entry_size)
   #         N = unpack(entry_format, record)[0]
   #         N = int(N)
    i = sum_n_vec(n)
    out = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids}
    return(out)

def myformula1(A,k,y,trA=None,trA2=None):
    std_y = (y-np.mean(y))/np.std(y)
    if (trA is None) and (trA2 is None):
        trA = np.sum(np.diag(A))
        trA2 = np.sum(np.multiply(A,A))
    n = A.shape[1]
    yay = std_y.T@A@std_y
    yty = std_y@std_y
    if (k==0):
        denominator = trA2 - 2*trA + n
        nominator = n - trA + yay - yty
        h2 = nominator/denominator
    else:
        pc = PCA_selected.values[:,:k]
        b = np.diag(pc.T@A@pc)-1
        c = (std_y@pc)**2 - 1
        denominator = trA2 - 2*trA + n - np.sum(b**2)
        nominator = n - trA + yay - yty - np.sum(b*c)
        h2 = nominator/denominator
    return(h2)

def myformula2(A,k,y,trA=None,trA2=None):
    y = y - np.mean(y)
    if (trA is None) and (trA2 is None):
        trA = np.sum(np.diag(A))
        trA2 = np.sum(np.multiply(A,A))
    n = A.shape[1]
    yay = y.T@A@y
    yty = np.dot(y,y)
    if (k==0):
        sigg = n*yay - trA*yty
        sige = trA2*yty - trA*yay
    else:
        pc = PCA_selected.values[:,:k]
        s = np.diag(pc.T@A@pc) #pciApci
        t = (y@pc)**2 #ypcipciy
        a11 = trA2 - np.sum(s**2)
        a12 = trA - np.sum(s)
        b1 = yay - np.sum(s*t)
        b2 = yty - np.sum(t)
        sigg = (n-k)*b1 - a12*b2
        sige = a11*b2 - a12*b1
    h2 = sigg/(sigg+sige)
    return(h2)

def regout(y,k):
    n = len(y)
    if (k==0):
        res = y
    else:
        pc = PCA_selected.values[:,:k]
        X = np.column_stack((np.ones(n),pc))
        XTX_inv = np.linalg.inv(np.dot(X.T,X))
        XTY = np.dot(X.T,y)
        beta = np.dot(XTX_inv,XTY)
        res = y - np.dot(X,beta)
    return(res)

G = ReadGRMBin(prefix)	
x = G['diag'].astype('float64')
n = x.size
#y = sp.spdiags(x, 0, x.size, x.size)
#GRM = csr_matrix(y)
idx = np.tril_indices(n,-1,n)
id_diag = np.diag_indices(n)
GRM_array_nona = np.zeros((n,n))
GRM_array_nona[idx] = G['off']
GRM_array_nona[id_diag] = x
GRM_array_nona.T[idx] = GRM_array_nona[idx]


trace_A = np.sum(x)
#trace_A2 = np.sum(np.multiply(GRM_array_nona,GRM_array_nona))
trace_A2 = 2*np.sum(G['off'].astype('float64')**2) + np.sum(x**2)
del(G)

PCA_selected = covariates.iloc[:,2:]

n_phen_nona = GRM_array_nona.shape[0]

y = phenotypes.values[:,2]
#f = open('stdy_'+output,'a')
#f.write(str(np.std(y))+'\n')
#f.close()
#std_y = (y-np.mean(y))/np.std(y)
#res_y = regout(std_y,n_covar)
#f = open('stdry_'+output,'a')
#f.write(str(np.std(res_y))+'\n')
#f.close()
res_y = regout(y,n_covar)
##Standardized AdjHE
h2_1 = myformula1(GRM_array_nona,n_covar,res_y,trA=trace_A,trA2=trace_A2)
end_time = timeit.default_timer()

f = open(output,'a')
f.write(str(h2_1)+'\n')
f.close()

#f = open('time'+output,'a')
#f.write(str(end_time-start_time)+'\n')
#f.close()
#


