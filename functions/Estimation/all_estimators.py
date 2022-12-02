#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 03:22:59 2022

@author: christian
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import itertools
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from functions.Data_input.load_data import load_everything
from functions.Estimation.AdjHE_estimator import load_n_AdjHE
from functions.Estimation.AdjHE_estimator import load_n_MOM
from functions.Estimation.PredLMM_estimator import load_n_PredLMM
from functions.Estimation.Estimate_helpers import create_formula
from functions.Estimation.GCTA_wrapper import GCTA


#%%
    

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
        result=  GCTA(df, covars, nnpc, mp, GRM, silent=False)
        
    if not silent :
        print(result["h2"])
    
    return pd.DataFrame(result, index = [0])

   
#%%
class Basu_estimation() :
    def __init__(self, prefix = None, pheno_file=None, cov_file=None, PC_file=None, k=0, ids = None, Simulation = False) :
        if prefix == None :
            print("Enter preloaded values...")
            self.df = None
            self.GRM = None
            self.phenotypes = "Y"
            self.simulation=True
        
        else : 
            print("Loading data...")
            self.df, self.GRM, self.phenotypes = load_everything(prefix, pheno_file, cov_file, PC_file, k, ids)
            self.simulation = False
    
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
            
    def estimate(self, npc, Method = None, RV = None, Naive = False, covars= False) : 
        print("Estimating")
        print(Method)
        if RV != None :
            print("RV: " + RV)
        # create empty list to store heritability estimates
        results = pd.DataFrame()
        
        
        # project GRM onto pc space
        # pcs = pd.DataFrame(PCA(n_components=20).fit_transform(self.GRM))
        # pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
        # # add the pcs to the dataframe
        # self.df = pd.concat([self.df, pcs], axis=1)


        # Forcing type to be integer for a little easier use
        if npc == None :
            npc = [0]
        
        if not covars :
            cov_combos = [[]]

        # Loop over each set of covariate combos
        for covs in self.cov_combos :
            # For each set of covariates recalculate the projection matrix
            
            # loop over all combinations of pcs and phenotypes
            for mp, nnpc in itertools.product(self.mpheno, npc):
                if not Naive :
                    r = load_n_estimate(
                        df=self.df, covars=covs, nnpc=nnpc, mp=mp, GRM= self.GRM, std= False, Method = Method, RV = RV)
                else :
                    # Empty results list
                    sub_results = pd.DataFrame({"Estimate": [],
                                           "Size": []})

                    # loop over  all sites
                    for site in np.unique(self.df[RV]):
                        # Grab the portion that lies within a given site
                        sub = self.df.loc[self.df.abcd_site == site, :]
                        # Get size
                        sub_n = sub.shape[0]
                        # Estimate just on the supsample
                        sub_result = load_n_estimate(df=sub, covars=[],  nnpc=0, mp= mp, GRM=self.GRM, std=False, Method= Method, RV=None,
                                                 silent=True)
                        sub_result = pd.DataFrame({"Estimate": [sub_result["h2"][0]],
                                               "Size": [sub_n]})
                        # Add to the list of estimates
                        sub_results = sub_results.append(sub_result, ignore_index=True)

                    # Pool the estimates
                    sub_results["nh2"] = (
                        sub_results["Size"] * sub_results["Estimate"]) / self.GRM.shape[0]
                    r = np.sum(sub_results["nh2"])
                    r = {"h2" : r, "ss": 0, "se" : 0, "var(sg)" : 0}
                    r = pd.DataFrame(r, index= [0])

                    
                    
                # Store each result
                results = pd.concat([results, r], ignore_index = True)


                
        self.results = results
        return self.results
    
    def pop_clusts(self, npc=2, groups = None):
        print("Generating PCA cluster visualization...")
        # Checked genotypes with coming from separate clusters
        trans = PCA(n_components=npc).fit_transform(self.GRM)

        if groups == None :
            sns.scatterplot(x=trans[:, 0], y=trans[:, 1])

        else: 
            sns.scatterplot(x=trans[:, 0], y=trans[:, 1],
                            hue=groups)

    def GRM_vis(self, sort_by=None, location=None, npc=0, plot_decomp = False):
        print("Generating GRM visualization...")
        if sort_by == None:
            df = self.df
        else:
            df = self.df.sort_values(sort_by)

        # Arrange the GRM in the same manner
        G = self.GRM[df.index, :][:, df.index]

        if npc != 0:
            # subtract substructre from GRM
            pcs = np.matrix(PCA(n_components=npc).fit(G).components_.T)
            P = pcs * np.linalg.inv(pcs.T * pcs) * pcs.T
            G2 = (np.eye(self.GRM.shape[0]) - P) * G
        else :
            G2 = G
        
        if plot_decomp :
            #subplot(r,c) provide the no. of rows and columns
            fig, ax = plt.subplots(1,3) 
            
            # use the created array to output your multiple images. In this case I have stacked 4 images vertically
            ax[0].imshow(G)
            ax[1].imshow(G2)
            P = P * G
            ax[2].imshow(P)
            ax[0].axis('off')  
            ax[1].axis('off')             
            ax[2].axis('off')             

            ax[0].set_title('GRM')
            ax[1].set_title('Residual relatedness')
            ax[2].set_title('Ethnicity contrib')      
            
        else: 
            plt.imshow(G2)
            # major_ticks = np.unique(df.abcd_site, return_counts= True)[1].cumsum()
            # plt.xticks(major_ticks -1)
            # plt.yticks(np.flip(major_ticks))
            # plt.grid(color='k', linestyle='-', linewidth=0.5)
            plt.axis('off')
    
            plt.title("GRM")
            plt.colorbar()

        if location != None:
            plt.savefig('docs/Presentations/Images/GRM_{n}_{s}_{c}_{site}.png'.format(
                n=self.nsubjects, s=self.nsites, c=self.nclusts, site=self.site_comp),
                dpi=200)

    
