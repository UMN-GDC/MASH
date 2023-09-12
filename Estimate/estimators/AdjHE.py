#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AdjHE estimator guts
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-05-26
"""

##############################################################
# These helpers stores all of the functions necessary for the 
# guts of the AdjHE estimator
##############################################################

import logging
import numpy as np
from scipy.linalg import block_diag
import pandas as pd



def AdjHE(A, df, mp, random_groups = None, npc=0, std=False):
    """
    Fucntion for generating heritability estimates from an Adjusted HE standpoint in closed form.
    Parameters
    ----------
    A : numpy array  
        n x n array containing the GRM values.
    df : pandas dataframe
        A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
    mp : string
        name of the phenotype to be estimated on.
    random_groups : string
        A string specifying what variable to use as a random variable. Default = None leads to no random variable
    npc : int, optional
        number of prinicpal components to adjust for. The default is 0.
    std : bool, optional
        Specify whether or not to standaridize the variables before estimation. The default is False.
    Returns
    -------
    tuple(ndarray, scalar)
        variance parameter estimates.
        standard error estimate
    """
    y = np.array(df[mp])
    y = (y-np.mean(y))/np.std(y)
    A = np.array(A)

    if npc ==0 :
        Sjs = 0
        Tjs =0
    else :
        # Grab the PCs
        PC_cols = [ col.startswith("pc")   for col in df ]
        # Select specified number of Pcs
        PC = np.array(df.iloc[:,PC_cols])[:, :npc]
        # Stack the PCs for easier matrix multiplications
        PCs = np.reshape(PC.T, newshape=(PC.shape[1], PC.shape[0], 1), order = "C")
        # transpose the PCs
        PCsT = np.reshape(PC.T, newshape=(PC.shape[1], 1, PC.shape[0]), order = "C")
    
    
        # Calculate outer products of eigenvectors
        PPt = np.matmul(PCs, PCsT)
        logging.debug(str(y.shape))
        logging.debug(str(PPt.shape))
        # Calculate tj's
        Tjs = np.matmul(np.matmul(y, PPt), y.T)
        # potentially could be Tjs = np.dot(y, PC) **2
        # Calculated Sj's
        Sjs = np.matmul(np.matmul(PCsT, A), PCs).flatten()
        
    # Compute elements of regression matrix
    n = A.shape[0]
    
    # Could square A faster with following A^2 = VD^2 V', but need to pass in eigenvalues D
    trA2 = np.trace(np.linalg.matrix_power(A, 2))
    trA = np.trace(A)
    
    topleft = trA2 - np.sum(Sjs ** 2)
    offdiag= trA - np.sum(Sjs)
    bottomright = n- npc
    
    if (random_groups is None) or (random_groups == "None") :
        # Solve the regression problem
        XXinv = np.linalg.inv(np.matrix([[topleft, offdiag],
                                         [offdiag, bottomright]]))
        yay = np.matmul(y.T, np.matmul(A, y))
        yty = np.inner(y,y)
        Ycol= np.array([yay - np.sum(Tjs * Sjs),
                        yty - np.sum(Tjs)])
        sigmas = abs(np.array(np.matmul(XXinv, Ycol)).flatten())
        
        h2 = sigmas[0] / (sigmas[0] + sigmas[1])
        # var_h2 = 2 * (sigmas[1]**2 * trA2 - 2*sigmas[0]*sigmas[1] * trA  + sigmas[0] **2 * n) / (sigmas[0] + sigmas[1])

        
    else :
        df = df.reset_index().drop("index", axis = 1)
        # Shuffle the A matrix so it matches the new order
        df = df.sort_values(random_groups).dropna(subset= [random_groups])
        A = np.matrix(A[df.index,:][:,df.index])
        n = A.shape[0]
        
        df["Intercept"] = 1
        # X_sites= df[rv]

        proj_cols = []
        # Grab the varibles for projection
        for col in df :
            if col.startswith("pc") :
                proj_cols.append(col)
        proj_cols += [random_groups]

        
        # Get dummies for categoricals if they exist
        X = np.matrix(pd.get_dummies(df[proj_cols]))

        # X = np.matrix(pd.get_dummies(df.drop(["FID", "IID", rv], axis = 1),  drop_first = True))
        y = np.matrix(df[mp])

        # Create S similarity matrix 
        site, sizes= np.unique(df[random_groups], return_counts = True)
        # Construct the block diagonal
        diags = [np.ones((size,size)) for size in sizes]
        S = np.matrix(block_diag(*diags))
        # Standardize S
        # S = (S - S.mean(axis = 1))/ S.std(axis = 1)

        # diags = [np.ones((size,size))* size for size in sizes]
        # S2 = np.matrix(block_diag(*diags) )
        
        # Construct the orthogonal projection matrix Q utilizing QR decomposition
        q, r = np.linalg.qr(X)
        Q = np.identity(n) - X.dot(np.linalg.inv(r).dot(q.T))
        Q = np.matrix(Q)
            
        # Compute elements of 3x3 matrix
        QSQ = Q * S * Q
        
        # find necessary traces
        trA2 = np.trace(A ** 2)
        trQSQA = np.trace(QSQ * A)
        trA = np.trace(A)
        trQSQ = np.trace(QSQ)
        trQSQQSQ = np.trace(QSQ * QSQ)

        # Find inverse 
        XtXm1 = np.linalg.inv(np.matrix([[trA2, trQSQA, trA],
                     [trQSQA, trQSQQSQ, trQSQ],
                     [trA, trQSQ, n]]))

        youter = np.matrix(np.outer(y,y))
        trAY = np.trace(A * youter)
        trQSQY = np.trace(QSQ * youter)
        trYout =  np.trace(youter)
        # Possible that the y's will need to account for prinicpal componetns in future real data cases
        sigmas = abs(XtXm1.dot(np.matrix([[trAY], [trQSQY], [trYout]])))
        
         
        h2 = sigmas[0,0] / (sigmas[0,0]  + sigmas[2,0])
        # var_h2 = 2 * (sigmas[2,0]**2 * trA2 - 2*sigmas[0,0]*sigmas[2,0] * trA  + sigmas[0,0] **2 * n) / (sigmas[0,0] + sigmas[2,0])

        
    # Calculate variance of estimate
    # var_h2 = 2/ ( trA2 - 2*trA + n - np.sum(Sjs**2))
    trA2 = np.trace(A ** 2)
    trA = np.trace(A)
    var_h2 = 2 / (trA2 - trA**2)
    if var_h2 < 0:
        logging.warning("Variance estimate is negative setting as absolute value")
        var_h2 = abs(var_h2)

    results = {"h2" : h2, "var(h2)" : var_h2}
    
    if results["h2"] < 0 :
        results["h2"] = 0
    elif results["h2"] >1 :
        results["h2"] = 1
    
    return results



