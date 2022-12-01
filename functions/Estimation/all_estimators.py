#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 03:22:59 2022

@author: christian
"""

import pandas as pd
import statsmodels.formula.api as smf
import itertools
from functions.Data_input.load_data import load_everything
from functions.Estimation.AdjHE_estimator import load_n_AdjHE
from functions.Estimation.AdjHE_estimator import load_n_MOM
from functions.Estimation.PredLMM_estimator import load_n_PredLMM
from functions.Estimation.Estimate_helpers import create_formula
from functions.Estimation.GCTA_wrapper import GCTA


#%%

class Basu_estimation() :
    def __init__(self, prefix, pheno_file, cov_file=None, PC_file=None, k=0, ids = None):
        # load data
        self.df, self.GRM, self.phenotypes = load_everything(prefix, pheno_file, cov_file, PC_file, k, ids)
    
    def looping(self, covars, npc, mpheno, loop_covars = False) :
        # Create list of covariate sets to regress over
        if covars != None :
            # Create the sets of covarates over which we can loop
            # This will return a list of lists of covariate names to regress on
            cov_combos = [covars[0:idx+1] for idx, c in enumerate(covars)]
            # If we don't want to loop, just grab the last item of the generated list assuming the user wants all of those variables included 
            if not loop_covars : 
                cov_combos = [cov_combos[-1]]
        else :
            cov_combos = [[]]
            
        self.cov_combos = cov_combos

        if mpheno == "all" :
            self.mpheno = self.phenotypes
        else :
            # make them lowercase
            self.mpheno =  mpheno
            
    def estimate(self, npc, Method = None, RV = None) : 
        # create empty list to store heritability estimates
        results = pd.DataFrame()

        # Forcing type to be integer for a little easier use
        if npc == None :
            npc = [0]

        # Loop over each set of covariate combos
        for covs in self.cov_combos :
            # For each set of covariates recalculate the projection matrix
            
            # loop over all combinations of pcs and phenotypes
            for mp, nnpc in itertools.product(self.mpheno, npc):
                r = load_n_estimate(
                    df=self.df, covars=covs, nnpc=nnpc, mp=mp, GRM= self.GRM, std= False, Method = Method, RV = RV)
                results = pd.concat([results, r], ignore_index = True)
                
        self.results
        return self.results


            
            
            
                
    

def load_n_estimate(df, covars, nnpc, mp, GRM, std = False, Method = "AdjHE", RV = None, silent=False):
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
        
    
        ids = df[["FID", "IID"]]
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
        nonmissing = ids[ids.IID.isin(temp.IID)].index
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
        result=  GCTA(df, covars, nnpc, mp, GRM, std = False, Method = "AdjHE", RV = None, silent=False)
        
    if not silent :
        print(result["h2"])
    
    return pd.DataFrame(result, index = [0])

   
    
