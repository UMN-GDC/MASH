#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 08:55:43 2022

@author: christian
"""

import numpy as np
#%%
                          
def columnwise_outter(mat) :
    """
    Gvien a (n x p) matrix of eigenvector loadings for individuals with p eigen loading vectors as the columns and n observations
    calculate and return all elementwise outer products stacked.

    Parameters
    ----------
    mat : numpy matrix
        (n x p) matrix with p columns for each eigenvetor loading and n observations.

    Returns
    -------
    St : np.array
        (p x n x n) stacked array of eigenvector loading outer products.

    """
    St = np.zeros((mat.shape[1], mat.shape[0], mat.shape[0]))

    for i in range(mat.shape[1]) :
        St[i,:,:] = np.outer(mat[:,i], mat[:,i])
    return St



def create_formula(nnpc, covars, mp, RV = None):
    """
    Creates a formula for the mean to residualize against. In other words describes the projection space and by extension the residual space.

    Parameters
    ----------
    nnpc : int
        number of PC's included in model.
    covars : list of integers
        the covariates to be included in the projection.
    mp : int
        phenotype to be projected.

    Returns
    -------
    (form, cols) where form is the formual exceptabule by smf.ols procedure and cols specifies the columns necessary for the projection
    (PCs, covariates, and phenotype).

    """
    # Get indices for ID variables
    id_cols = ["fid", "iid"] 
    # Get the full range of pc columns
    pc_cols = ["pc_" + str(p) for p in range(1, nnpc +1)]
    # Create formula string
    
    if len(covars) != 0:
        RHS = " + ".join(covars)
    else :
        RHS = "1"


    form = mp + "~ " +  RHS
    # columns
    cols = id_cols + [mp] + covars + pc_cols
    if RV != None :
        cols = cols + [RV]

    # return the formula and columns
    return(form, cols)
    
    
