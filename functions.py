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


def myformula1(A,y,trA=None,trA2=None, npc=0):
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


def myformula2(A,y,trA=None,trA2=None, npc =0 ):
#    y = y - np.mean(y)
    if (trA is None) and (trA2 is None):
        trA2 = np.sum(np.multiply(A,A))
        trA = np.sum(np.diag(A))

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

# only if covaraites or PC's nonnull
def regout(X, y, covar, PC):
    if (covar!="NULL") and (PC != "NULL"):
        X = cov_selected
        XTX_inv = np.linalg.inv(np.dot(X.T,X))
        XTY = np.dot(X.T,y)
        beta = np.dot(XTX_inv,XTY)
        res = y - np.dot(X,beta)
        return(res)
    else: 
        return(y)


