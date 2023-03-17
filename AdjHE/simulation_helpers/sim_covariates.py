#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:56:41 2023

@author: christian
"""
import numpy as np
import pandas as pd
import scipy

def sim_covariates(rng, nclusts, df, cov_effect=True, ortho_cov = False): 
    """
    Generate the covariates as well as associated effects for each subject

    Parameters
    ----------
    rng : numpy random number generator object
        numpy random number genrator object.
    nclusts : int
        number of genetic clusters in the dataset.
    df : pandas datafraem
        dataframe containing information on each subject for the simulation.
    cov_effect : bool, optional
        Should the covariates cause a fixed effect or no effect at all?. The default is True.
    ortho_cov : bool, optional
        should the covariates be orthogonal to principal components. This is generally overkill for high dimensional data where
        n>> p. The default is False.

    Returns
    -------
    Xc : numpy array
        the observed varaible for each subject.
    Covar_contrib : numpy array
        the contribution to the phenotype from the covaraites.

    """
    # To generate orthogonal, find the Null space of the transpose of the principal components and site dummy matrix
    # Then, construct a random vector from the null space.
    nsubjects = df.shape[0]
    if ortho_cov :
        pcs = ["pc_" + str(i + 1) for i in range(nclusts)]
        cov_arr= pd.get_dummies(df["abcd_site"].astype(str)).to_numpy()
        cov_arr = np.concatenate((cov_arr, df[pcs].to_numpy()), axis =1 )
        ns = scipy.linalg.null_space(cov_arr.T)
        Xc = np.dot(ns, np.random.rand(ns.shape[1],1))

    else: 
        Xc = rng.uniform(0,1, nsubjects)
    
    
    if cov_effect :
        # Scale it to be on same size as genetic contribution
        Covar_contrib = Xc * np.mean(df["Gene_contrib"])
    else :
        df["Covar_contrib"] = 0 
        
    return Xc, Covar_contrib

