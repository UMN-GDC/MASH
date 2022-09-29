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
    data = data.sort_values(rv)
    A = A[data.index,:][:,data.index]
    n = A.shape[0]
    
    X = data.drop(mp, axis =1)
    y = data[mp]

    # Create S similarity matrix 
    site, sizes= np.unique(data[rv], return_counts = True)
    #%% Construct the block diagonal
    diags = [np.ones((size,size)) for size in sizes]
    S = block_diag(*diags)
    diags = [np.ones((size,size)) for size in sizes]

    
    # Construct the orthogonal projection matrix Q utilizing QR decomposition
    q, r = np.linalg.qr(X)
    Q = np.identity(n) - X.dot(np.linalg.inv(r).dot(q.T))

        
    # Compute elements of 3x3 matrix
    QS = Q.dot(S)
    QSQ = QS.dot(Q)
    QA = Q.dot(A)
    
    # find necessary traces
    trA2 = np.trace(np.linalg.matrix_power(A,2))
    trQSQA = np.trace(QS.dot(QA))
    trA = np.trace(A)
    trQS2Q= np.trace(np.linalg.matrix_power(QS, 2))
    trQSQ = np.trace(QSQ)
    
    
    # Find solution 
    XtXm1 = np.linalg.inv(np.matrix([[trA2, trQSQA, trA],
                     [trQSQA, trQS2Q, trQSQ],
                     [trA, trQSQ, n]]))
    sigmas = np.dot(XtXm1, np.concatenate([A, QSQ, np.identity(n)], axis = 0) )* y 
    
    return(sigmas)
    
    
   
   
    
