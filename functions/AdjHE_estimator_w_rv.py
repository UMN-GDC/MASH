#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:55:43 2022

@author: christian
"""

import numpy as np
from scipy.linalg import block_diag
import pandas as pd
#%%


def AdjHE_rv_estimator(A,data, mp, rv, npc=0, std=False) :
    """
    Estimate the heritability of the presence of an additional random effect.

    Parameters
    ----------
    A : numpy array  
        n x n array containing the GRM values.
    data : pandas dataframe
        A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
    mp : int
        1 based index specifying the phenotype to be estimated on.
    rv : string
        specifying the name of the column to be used as a random variable.
    npc : int, optional
        number of prinicpal components to adjust for. The default is 0.
    std : bool, optional
        Specify whether or not to standaridize the variables before estimation. The default is False.

    Returns
    -------
    tuple(scalar, scalar)
        h2 - heritability estimate.
        standard error estimate
    """
    # Reorder df by the random variable
    # then reorder the GRM to match
    print("AdjHE + random effect")
    data = data.reset_index().drop("index", axis = 1)
    data = data.sort_values(rv).dropna(subset= [rv])
    A = np.matrix(A[data.index,:][:,data.index])
    n = A.shape[0]
    
    data["Intercept"] = 1
    X = np.matrix(data.drop([mp, "fid", "iid", rv], axis =1))
    y = np.matrix(data[mp])

    # Create S similarity matrix 
    site, sizes= np.unique(data[rv], return_counts = True)
    #%% Construct the block diagonal
    diags = [np.ones((size,size)) for size in sizes]
    S = np.matrix(block_diag(*diags))
    diags = [np.ones((size,size))* size for size in sizes]
    S2 = np.matrix(block_diag(*diags) )
    
    # Construct the orthogonal projection matrix Q utilizing QR decomposition
    q, r = np.linalg.qr(X)
    Q = np.identity(n) - X.dot(np.linalg.inv(r).dot(q.T))
    Q = np.matrix(Q)
        
    # Compute elements of 3x3 matrix
    QS = Q * S
    QSQ = QS * Q
    
    # find necessary traces
    trA2 = np.trace(A ** 2)
    trQSQA = np.trace(QSQ * A)
    trA = np.trace(A)
    trQSQ = np.trace(QSQ)
    trQS2Q = np.trace(Q * S2 * Q)
    
    youter = np.matrix(np.outer(y,y))
    trAY = np.trace(A * youter)
    trQSQY = np.trace(QSQ * youter)
    trYout = np.linalg.norm(y) * n
    # Find solution 
    XtXm1 = np.linalg.inv(np.matrix([[trA2, trQSQA, trA],
                     [trQSQA, trQS2Q, trQSQ],
                     [trA, trQSQ, n]]))
    # Possible that the y's will need to account for prinicpal componetns in future real data cases
    sigmas = XtXm1.dot(np.matrix([[trAY], [trQSQY], [trYout]]))
    
    # return heritability estimate
    return (sigmas[0] / sum(sigmas))[0,0], 0
    
    
   
   
    
