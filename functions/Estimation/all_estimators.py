#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 03:22:59 2022

@author: christian
"""

import pandas as pd
import statsmodels.formula.api as smf
from functions.Estimation.AdjHE_estimator import load_n_AdjHE
from functions.Estimation.AdjHE_estimator import load_n_MOM
from functions.Estimation.PredLMM_estimator import load_n_PredLMM
from functions.Estimation.Estimate_helpers import create_formula
from functions.Estimation.GCTA_wrapper import GCTA


#%%

def load_n_estimate(df, covars, nnpc, mp, GRM, std = False, Method = "AdjHE", RV = None, silent=False, args = None):
    """
    Estimates heritability, but solves a full OLS problem making it slower than the closed form solution. Takes 
    a dataframe, selects only the necessary columns (so that when we do complete cases it doesnt exclude too many samples)
    residualizes the phenotype, then documents the heritability, standard error and some computer usage metrics.

    Parameters
    ----------
    df : pandas dataframe
        dataframe contianing phenotype, covariates, an prinicpal components.
    covars : list of int
        list of integers specifying which covariates to include in the resiudalization.
    nnpc : int
        number of pcs to include.
    mp : int
        which phenotype to estiamte on.
    GRM : np array
        the GRM with missingness removed.
    std : bool, optional
        specifying whether standarization happens before heritability estimation. The default is False.
    Method: str
        specify which method of estimation to use AdjHE, PredLMM, or MOM
        Default is AdjHE
    RV : string, optional
        Varible to control for as a random effect, if applicable

    Returns
    -------
    pandas dataframe containing:
        - heritability estimate
        - standard error the estimate
        - the phenotype
        - the number of pcs included
        - The covarites included 
        - time for analysis
        - maximum memory usage
    """
    if not silent :
        print(Method + "Estimation...")

    # Remove missingness for in-house estimators
    if Method != "GCTA" :
        
    
        ids = df[["fid", "iid"]]
        # seed empty result vector
        # result.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]
        # create the regression formula and columns for seelcting temporary
        form, cols  = create_formula(nnpc, covars, mp, RV)
        # save a temporary dataframe
        temp = df[cols].dropna()
        # Save residuals of selected phenotype after regressing out PCs and covars
        temp[mp] = smf.ols(formula = form, data = temp, missing = 'drop').fit().resid
        # Potentially could use this to control for random effects
        # smf.mixedlm(formula= form, data = temp, groups=temp["scan_site"])
        # keep portion of GRM without missingess for the phenotypes or covariates
        nonmissing = ids[ids.iid.isin(temp.iid)].index
        GRM_nonmissing = GRM[nonmissing,:][:,nonmissing]
        print(temp.columns)

    # Select method of estimation
    if Method == "AdjHE": 
        result = load_n_AdjHE(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)

    elif Method == "MOM": 
        result = load_n_MOM(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)
    elif Method == "PredlMM" : 
        result = load_n_PredLMM(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)
    elif Method == "GCTA" :
        
        result = GCTA(grm = args["prefix"], pheno_file = args["pheno"], cov_file = args["covar"], PC_file = args["PC"], covars = covars, 
                      nnpc = nnpc, mp = mp)

    if not silent :
        print(result["h2"])
    
    return(pd.DataFrame(result))

   
    
