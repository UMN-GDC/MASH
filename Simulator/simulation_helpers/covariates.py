#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:56:41 2023

@author: christian
"""
import numpy as np
import pandas as pd
import scipy

def sim_covariates(rng, nclusts, df): 
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

    Xc = rng.uniform(0,1, nsubjects)
    
    Covar_contrib = Xc * np.mean(df["resid_eff"])
        
    return Xc, Covar_contrib

