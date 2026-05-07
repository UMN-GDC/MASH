#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 08:31:35 2022

@author: christian
"""
#import os
#os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
import numpy as np
import statsmodels.formula.api as smf
import pandas as pd

#%%
def first_level_reg(data, dep_var, fixed) :
    """
    Performs a first level regression of the fixed effects.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe containing the necessary fixed effects as well as the phenotype of interest.
    dep_var : string
        String specifying the column to use as a dependent variable.
    fixed : string, or list of strings
        The set of variables to treat as fixed effectors of the dependent variable.

    Returns
    -------
    vector of residuals from model fit.

    """
    # Create the first level formula
    ## If it's a string just make it the RHS
    if isinstance(fixed, str) :
        RHS = fixed
    ## If it's a list join them
    elif isinstance(fixed, list) :
        RHS = " + ".join(fixed)
    ## Otherwise just control for the intercept
    else: 
        RHS = "1"

    # Do the actual first level regression
    form = " ~ ".join([dep_var, RHS])
    model = smf.ols(form, data = data).fit()
    # print(model1.summary())
    
    return model.resid

def second_level_reg(resids, GRM, RV = None) : 
    """
    2nd level regression of the residuals of the first model vs the GRM and the second specified random variable

    Parameters
    ----------
    data : pandas dataframe
        Dataframe conatining dep endent varaible, and fixed effect variables
    model1 : statsmodels fit model object
        Fit model object from the first level regression.
    GRM : np.array 2D
        Genetic Relatedness Matrix nxn.
    RV : np.array 2D
        outer product between the second specified random variable.

    Returns
    -------
    model : statsmodels fit model object
        fit model from the second level regression.

    """
    # Get the number of observations
    n = GRM.shape[0] - 1
    # Create a dataframe for all of the variable and outer products, only keeping the lower half of the matrix, since they are 
    # symmetric
    df = pd.DataFrame({
        "yy" : np.outer(resids,  resids)[np.tril_indices(n)],
        "GRM" : np.array(GRM)[np.tril_indices(n)],
        "II" : np.identity(n)[np.tril_indices(n)]})
    # Allow option to not specify additional RV
    if RV != None :
        df[RV] = RV[np.tril_indices(n)]
        model = smf.ols("yy ~ GRM + RV + II-1 ", data = df).fit()
    else :
        model = smf.ols("yy ~ GRM + II-1 ", data = df).fit()
    return model


def two_level_regression(data, GRM, fixed, dep_var, RV = None) :
    # fit the first level model
    m1 = first_level_reg(data, dep_var, fixed)
    # pass the first model fit to the second level
    m2 = second_level_reg(m1, GRM, RV)
    return m2


