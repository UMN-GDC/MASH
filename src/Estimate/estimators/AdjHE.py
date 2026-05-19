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
    if std:
        y = (y - y.mean()) / y.std()
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
        PCs = np.reshape(PC.T, (PC.shape[1], PC.shape[0], 1), order="C")
        # transpose the PCs
        PCsT = np.reshape(PC.T, (PC.shape[1], 1, PC.shape[0]), order="C")
    
        # Calculate outer products of eigenvectors
        PPt = np.matmul(PCs, PCsT)
        # Calculate tj's
        Tjs = np.matmul(np.matmul(y, PPt), y.T)
        # Calculated Sj's
        Sjs = np.matmul(np.matmul(PCsT, A), PCs).flatten()
        
    # Compute elements of regression matrix
    n = A.shape[0]
    
    trA2 = np.trace(np.linalg.matrix_power(A, 2))
    trA = np.trace(A)
    
    topleft = trA2 - np.sum(Sjs ** 2)
    offdiag= trA - np.sum(Sjs)
    bottomright = n- npc
    
    if random_groups is None :
        # Solve the regression problem
        XXinv = np.linalg.inv(np.array([[topleft, offdiag],
                                         [offdiag, bottomright]]))
        yay = np.matmul(y.T, np.matmul(A, y))
        yty = np.inner(y,y)
        Ycol= np.array([yay - np.sum(Tjs * Sjs),
                        yty - np.sum(Tjs)])
        sigmas = np.array(XXinv @ Ycol).flatten()
        
        h2 = sigmas[0] / (sigmas[0] + sigmas[1])
        if sigmas[0] < 0:
            logging.warning(f"AdjHE: Negative genetic variance. sigmas=({sigmas[0]:.4e}, {sigmas[1]:.4e}), "
                           f"yay={yay:.4e}, yty={yty:.4e}, trA2={trA2:.4e}, trA={trA:.4e}, "
                           f"topleft={topleft:.4e}, offdiag={offdiag:.4e}, bottomright={bottomright:.4e}")
        ss = 0
        # var_h2 = 2 * (sigmas[1]**2 * trA2 - 2*sigmas[0]*sigmas[1] * trA  + sigmas[0] **2 * n) / (sigmas[0] + sigmas[1])

        
    else :
        # Shuffle the A matrix so it matches the new order
        df = df.sort_values(random_groups).dropna(subset= [random_groups])
        A = np.array(A[df.index,:][:,df.index])
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
        X = np.array(pd.get_dummies(df[proj_cols]))

        # X = np.matrix(pd.get_dummies(df.drop(["FID", "IID", rv], axis = 1),  drop_first = True))
        y = np.array(df[mp])
        if std:
            y = (y - y.mean()) / y.std()

        # Create S similarity matrix 
        site, sizes= np.unique(df[random_groups], return_counts = True)
        # Construct the block diagonal
        diags = [np.ones((size,size)) for size in sizes]
        S = np.array(block_diag(*diags))
        # Standardize S
        # S = (S - S.mean(axis = 1))/ S.std(axis = 1)

        # diags = [np.ones((size,size))* size for size in sizes]
        # S2 = np.matrix(block_diag(*diags) )
        
        # Construct the orthogonal projection matrix Q utilizing QR decomposition
        q, r = np.linalg.qr(X)
        Q = np.identity(n) - X @ np.linalg.inv(r) @ q.T
        Q = np.array(Q)
            
        # Compute elements of 3x3 matrix
        QSQ = Q @ S @ Q
        
        # find necessary traces
        trA2 = np.trace(A @ A)
        trQSQA = np.trace(QSQ @ A)
        trA = np.trace(A)
        trQSQ = np.trace(QSQ)
        trQSQQSQ = np.trace(QSQ @ QSQ)

        # Find inverse 
        XtXm1 = np.linalg.inv(np.array([[trA2, trQSQA, trA],
                     [trQSQA, trQSQQSQ, trQSQ],
                     [trA, trQSQ, n]]))

        youter = np.outer(y,y)
        trAY = np.trace(A @ youter)
        trQSQY = np.trace(QSQ @ youter)
        trYout =  np.trace(youter)
        # Possible that the y's will need to account for prinicpal componetns in future real data cases
        sigmas = XtXm1 @ np.array([[trAY], [trQSQY], [trYout]])
        
        if sigmas[0,0] < 0:
            logging.warning(f"AdjHE (random_groups): Negative genetic variance. "
                           f"sigmas=({sigmas[0,0]:.4e}, {sigmas[1,0]:.4e}, {sigmas[2,0]:.4e})")
        
        ss = sigmas[1,0]
        
        h2 = sigmas[0,0] / (sigmas[0,0] + sigmas[2,0])
        # var_h2 = 2 * (sigmas[2,0]**2 * trA2 - 2*sigmas[0,0]*sigmas[2,0] * trA  + sigmas[0,0] **2 * n) / (sigmas[0,0] + sigmas[2,0])

        
    # Calculate variance of estimate
    trA2 = np.trace(np.linalg.matrix_power(A, 2))
    trA = np.trace(A)
    var_h2 = 2 * n / (n * trA2 - trA**2)
    if var_h2 <= 0:
        logging.warning(f"AdjHE: Variance estimate is negative ({var_h2:.4e}), setting as absolute value")
        var_h2 = abs(var_h2)

    h2_raw = h2
    results = {"h2" : h2, "var(h2)" : var_h2, "h2_raw" : h2_raw}
    
    if results["h2"] < 0 :
        logging.warning(f"AdjHE: h2 clamped from {results['h2']:.4e} to 0. "
                       f"This indicates heritability near zero or a poor model fit")
        results["h2"] = 0
    elif results["h2"] > 1 :
        logging.warning(f"AdjHE: h2 clamped from {results['h2']:.4e} to 1")
        results["h2"] = 1
    
    return results

#%%



# def AdjHE_rv_estimator_homo(A,df, mp, rv, npc=0, std=False) :
#     """
#     Estimate the heritability of the presence of an additional random effect.

#     Parameters
#     ----------
#     A : numpy array  
#         n x n array containing the GRM values.
#     df : pandas dataframe
#         A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
#     mp : int
#         1 based index specifying the phenotype to be estimated on.
#     rv : string
#         specifying the name of the column to be used as a random variable.
#     npc : int, optional
#         number of prinicpal components to adjust for. The default is 0.
#     std : bool, optional
#         Specify whether or not to standaridize the variables before estimation. The default is False.

#     Returns
#     -------
#     tuple(scalar, scalar)
#         h2 - heritability estimate.
#         standard error estimate
#     """
#     # Reorder df by the random variable
#     # then reorder the GRM to match
#     df = df.reset_index().drop("index", axis = 1)
#     df = df.sort_values(rv).dropna(subset= [rv])
#     A = np.matrix(A[df.index,:][:,df.index])
#     n = A.shape[0]
    
#     df["Intercept"] = 1
#     X_sites= df[rv]

#     # Get dummies for categoricals if they exist
#     X = np.matrix(pd.get_dummies(df.drop(["FID", "IID", rv], axis = 1),  drop_first = True))
#     y = np.matrix(df[mp])

#     # Create S similarity matrix 
#     site, sizes= np.unique(df[rv], return_counts = True)
#     # Construct the block diagonal
#     diags = [np.ones((size,size)) for size in sizes]
#     S = np.matrix(block_diag(*diags))
#     # Standardize S
#     # S = (S - S.mean(axis = 1))/ S.std(axis = 1)

#     # diags = [np.ones((size,size))* size for size in sizes]
#     # S2 = np.matrix(block_diag(*diags) )
    
#     # Construct the orthogonal projection matrix Q utilizing QR decomposition
#     q, r = np.linalg.qr(X)
#     Q = np.identity(n) - X.dot(np.linalg.inv(r).dot(q.T))
#     Q = np.matrix(Q)
        
#     # Compute elements of 3x3 matrix
#     QSQ = Q * S * Q
    
#     # find necessary traces
#     trA2 = np.trace(A ** 2)
#     trQSQA = np.trace(QSQ * A)
#     trQSQQSQ = np.trace(QSQ * QSQ)


