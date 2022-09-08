#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:55:43 2022

@author: christian
"""

import numpy as np
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
    rv : int
        1 based index specifying the variable to be used as a random variable.
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
    y = data[mp]
    # Extract random variable
    RV = data[rv]
    # Create S similarity matrix 
    S = np.outer(RV, RV)
        
    # Compute elements of 3x3 matrix
    trA2 = np.tr(A**2)
    trSA = np.tr(S*A)
    trA = np.tr(A)
    n = A.shape[0]
    
    
    # Find solution 
    XtX = np.matrix([[trA2, trSA, trA],
                     [trSA, n**2, n],
                     [trA, n, n]])
    sigmas = np.linalg.inv(XtX) * np.concatenate([A, S, np.identity(n)], axis = 0) * y 
    
    return(sigmas)
    
    
   
   
    
   
#%%    

A = GRM_array_nona
mp = 1
rv = 4
data = df
npc= 0
std = False

#%%

