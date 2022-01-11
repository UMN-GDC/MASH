import pandas as pd
import numpy as np
import sys
import os
import psutil
import scipy.sparse as sp
import scipy.sparse.linalg
import inspect
from scipy.sparse import csr_matrix, rand
from struct import unpack, calcsize
from numpy.linalg import inv
from numpy.linalg import multi_dot
import timeit
import logging
import resource
import argparse
from argparse import RawTextHelpFormatter
from functions import *



def sum_n_vec(n):
    s = [int(0)] * n
    for i in range(n):
        s[i] = int(((i + 1) * (i + 2) / 2) - 1)
    return(s)


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
    val = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids,'N':N}
    return(val)


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
    val = {'a':incr,'b':np.repeat(counts,counts)}
    return(val)




# Read data function that can load csv pheno and txt file types
def read_datas(file_path, GRM_ids) :
 if(file_path.split(".")[-1] == "csv"):
  dat = pd.read_csv(file_path)
 elif(file_path.split(".")[-1] == "phen"):
  dat = pd.read_table(file_path, sep = " " )
 elif(file_path.split(".")[-1] == "txt"):
  dat = pd.read_table(file_path, sep = " " )
 # remove the unintentional columns that sometimes happen with phenotype and csv filetypes
 dat = dat[dat.columns.drop(list(dat.filter(regex='Unnamed')))]
 # only keep intersections of the individuals in the GRM and in the covariates and phenotypes
 intersection_indiv = np.intersect1d(GRM_ids.iloc[:,1], dat.IID)
 dat = dat[dat.IID.isin(intersection_indiv)]
 # store it as an array for efficiency and for certain linear alg functions to work
 dat = dat.drop(["FID", "IID"], axis = 1)
# dat = np.asarray(dat.drop(["FID", "IID"], axis = 1))
 return(dat)

