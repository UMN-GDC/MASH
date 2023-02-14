#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 03:22:59 2022

@author: christian
"""
import timeit
import resource
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import itertools
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
from functions.Data_input.load_data import load_everything
from functions.Estimation.AdjHE_estimator import AdjHE_estimator #, load_n_MOM
from functions.Estimation.PredLMM_estimator import load_n_PredLMM
from functions.Estimation.GCTA_wrapper import gcta, GCTA
from functions.Estimation.combat import neuroCombat


# %%

def load_n_estimate(df, fixed_effects, nnpc, mp, GRM, std=False, Method="AdjHE", random_groups=None, silent=False, homo=True, gcta=gcta):
    """
    Estimates heritability, but solves a full OLS problem making it slower than the closed form solution. Takes 
    a dataframe, selects only the necessary columns (so that when we do complete cases it doesnt exclude too many samples)
    residualizes the phenotype, then documents the heritability, standard error and some computer usage metrics.

    Parameters
    ----------
    df : pandas dataframe
        dataframe contianing phenotype, covariates, an prinicpal components.
    fixed_effects : list of int
        list of integers specifying which covariates to include in the resiudalization. This should not include the random_groups
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
    random_groups : string, optional
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

    # Remove missingness for in-house estimators
    if Method != "GCTA":
        id_cols = ["FID", "IID"]
        ids = df[id_cols]

        # Change nnpc to a number
        if nnpc == None :
            nnpc =0
            
        pc_cols = ["pc_" + str(p) for p in range(1, nnpc +1)]

        if fixed_effects == None :
            fixed_effects = []
        

        
        # Create formula string
        if len(fixed_effects) != 0:
            RHS = " + ".join(fixed_effects)
        else :
            RHS = "1"
        
        # Take care of making random groups a list
        if random_groups == None :
            group_cols = []
        elif isinstance(random_groups, str) :
            group_cols = [random_groups]

        # columns
        cols = id_cols + [mp] + fixed_effects + pc_cols + group_cols
        
        # Make formula
        form = mp + "~ " +  RHS
        
        # save a temporary dataframe
        temp = df[cols].dropna()
        nonmissing = ids[ids.IID.isin(temp.IID)].index
        GRM_nonmissing = GRM[nonmissing, :][:, nonmissing]
        
    # Select method of estimation
    if Method == "AdjHE":
        # AdjHE projects away covariates to start
        temp[mp] = smf.ols(formula=form, data=temp, missing='drop').fit().resid
        result = AdjHE_estimator(A = GRM_nonmissing, df=temp, mp = mp, random_groups = random_groups, npc= nnpc, std=std)

    # MOM estimator is under construction
    # elif Method == "MOM":
    #     result = load_n_MOM(temp, covars, nnpc, mp,
    #                         GRM_nonmissing, std=False, RV=RV)
    elif Method == "PredlMM":
        result = load_n_PredLMM(temp, fixed_effects, nnpc, mp,
                                GRM_nonmissing, std=False, random_groups=random_groups)
        
    elif Method == "GCTA":
        result = GCTA(df, fixed_effects + [], nnpc, mp, GRM, gcta=gcta, silent=False)
        
    elif Method == "SWD":
        # SWD projects away sites then projects away covaraites
        temp[mp] = smf.ols(formula= mp + " ~ " + random_groups, data=temp, missing='drop').fit().resid
        temp[mp] = smf.ols(formula=form, data=temp, missing='drop').fit().resid        
        result = AdjHE_estimator(A = GRM_nonmissing, df = temp, mp = mp, random_groups = None, npc=nnpc, std=False)

        
    elif Method == "Combat":
        # Format data for Harmonization tool
        temp["Y1"] = df["Y1"]
        tempy = temp[["Y", "Y1"]].T

        # Harmonization step:
        data_combat = neuroCombat(dat=tempy,
                                  covars=df[["pc_" + str(i + 1) for i in range(nnpc)] + fixed_effects],
                                  batch_col=random_groups)["data"]

        # Grab the first column of the harmonized data
        temp[mp] = data_combat[0, :].T
        result = AdjHE_estimator(A = GRM_nonmissing, df = temp, mp = mp, random_groups = None, npc=nnpc, std=False)
    

    else:
        print("Not an accepted method")
        
    

    return pd.DataFrame(result, index=[0])


# %%
class Basu_estimation():
    def __init__(self, prefix=None, pheno_file=None, cov_file=None, PC_file=None, k=0, ids=None, Simulation=False):
        if prefix == None:
            print("Enter preloaded values...")
            self.df = None
            self.GRM = None
            self.phenotypes = "Y"
            self.simulation = True

        else:
            print("Loading data...")
            self.df, self.GRM, self.phenotypes = load_everything(
                prefix, pheno_file, cov_file, PC_file, k, ids)
            self.simulation = False

    def estimate(self, npc, mpheno="all", Method=None, random_groups =None, Naive=False, fixed_effects=None, homo=True, loop_covars=False):
        
        start_est = timeit.default_timer()
                
        # Create list of covariate sets to regress over
        if (fixed_effects == None) or (fixed_effects == []):
            fixed_combos = [[]]
        else:
            # Create the sets of covarates over which we can loop
            # This will return a list of lists of covariate names to regress on
            fixed_combos = [fixed_effects[0:idx+1] for idx, c in enumerate(fixed_effects)]
            # If we don't want to loop, just grab the last item of the generated list assuming the user wants all of those variables included
            if not loop_covars:
                fixed_combos = [fixed_combos[-1]]

        if mpheno == "all":
            self.mpheno = self.phenotypes
        else:
            # make them lowercase
            self.mpheno = mpheno

        print("Estimating")
        print(Method)
        if random_groups != None:
            print("RV: " + random_groups)
        # create empty list to store heritability estimates
        results = pd.DataFrame()

        # project GRM onto pc space
        # pcs = pd.DataFrame(PCA(n_components=20).fit_transform(self.GRM))
        # pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
        # # add the pcs to the dataframe
        # self.df = pd.concat([self.df, pcs], axis=1)

        # Forcing type to be integer for a little easier use
        if npc == None:
            npc = [0]

        # Loop over each set of covariate combos
        for covs in tqdm(fixed_combos, desc = "Covariate sets"):
            # For each set of covariates recalculate the projection matrix

            # loop over all combinations of pcs and phenotypes
            for mp, nnpc in tqdm(itertools.product(self.mpheno, npc), desc = "Phenotype, PC combination counter"):
                
                start_est = timeit.default_timer()
                try: 
                    C = "+".join(fixed_combos)
                except TypeError :
                    C = "None"
                    
                r = {"Pheno": mp, 
                          "PCs" : nnpc,
                          "Covariates" : C}

                if not Naive:
                    r = load_n_estimate(
                        df=self.df, fixed_effects=covs, nnpc=nnpc, mp=mp, GRM=self.GRM, std=False, Method=Method, random_groups=random_groups, homo=homo)

                else:
                    # Empty results list
                    sub_results = pd.DataFrame({"h2": [],
                                                "Size": []})

                    # loop over  all sites
                    for group in np.unique(self.df[random_groups]):
                        # Grab the portion that lies within a given site
                        sub_df = self.df.loc[self.df[random_groups] == group, :].reset_index(drop = True)
                        # Get size
                        sub_n = sub_df.shape[0]
                        sub_GRM = self.GRM[self.df[random_groups] == group,:][:,self.df[random_groups] == group]
                        # Find PC's individually for each site
                        if nnpc != 0:
                            pcs = pd.DataFrame(PCA(n_components=10).fit_transform(np.asarray(sub_GRM)))
                            pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
                            # Drop previous pcs
                            sub_df = sub_df.loc[:,~sub_df.columns.str.startswith('pc_')]
                        #Add site specific PC's
                            sub_df = pd.concat([sub_df, pcs], axis=1)



                        # Estimate just on the supsample
                        sub_result = load_n_estimate(df=sub_df, fixed_effects=[],  nnpc=nnpc, mp=mp, GRM=sub_GRM, std=False, Method=Method, random_groups=None,
                                                     silent=True, homo=homo)
                        sub_result = pd.DataFrame({"h2": [sub_result["h2"][0]],
                                                   "Size": [sub_n]})
                        # Add to the list of estimates
                        sub_results = sub_results.append(
                            sub_result, ignore_index=True)

                    # Pool the estimates
                    sub_results["nh2"] = (
                        sub_results["Size"] * sub_results["h2"]) / self.GRM.shape[0]
                    h2 = np.sum(sub_results["nh2"])
                    r["h2"] = h2
                    r["ss"] = np.nan
                    r["var(h2)"] =  np.var(sub_results["h2"])
                    r = pd.DataFrame(r, index=[0])

                # Store each result
                results["h2"] = r["h2"]
                results["ss"] = r["ss"]
                results["var(h2)"] = r["var(h2)"]
                results["Analysis time"] = timeit.default_timer() - start_est
                # Get memory for each step (in Mb) (This is a little sketchy)
                results["mem"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000

                results = pd.concat([results, r], ignore_index=True)

        self.results = results
        return results # , sub_results

    def pop_clusts(self, npc=2, groups=None):
        print("Generating PCA cluster visualization...")
        # Checked genotypes with coming from separate clusters
        pca = PCA(n_components=npc).fit(self.GRM)
        trans = pca.transform(self.GRM)

        sns.scatterplot(x=trans[:, 0], y=trans[:, 1],
                        hue=self.df[groups].astype(str))
        plt.show()

        plt.plot(pca.explained_variance_ratio_)
        plt.show()

    def GRM_vis(self, sort_by=None, location=None, npc=0, plot_decomp=False):
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
        else:
            G2 = G

        if plot_decomp:
            # subplot(r,c) provide the no. of rows and columns
            fig, ax = plt.subplots(1, 3)

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
