#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 03:22:59 2022

@author: christian
"""
import timeit
import logging
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import itertools
from sklearn.decomposition import PCA

from tqdm.auto import tqdm
from Estimate.data_input.load_data import load_everything
from Estimate.estimators.AdjHE import AdjHE #, load_n_MOM


def sort_pc_columns(pc_cols):
    """Sort PC column names numerically (pc1, pc2, ..., pc10, pc11)."""
    def get_pc_num(col):
        import re
        match = re.search(r'pc(\d+)', col.lower())
        return int(match.group(1)) if match else float('inf')
    return sorted(pc_cols, key=get_pc_num)
from Estimate.estimators.PredLMM import load_n_PredLMM
from Estimate.estimators.GCTA_wrapper import gcta, GCTA
from Estimate.estimators.combat import neuroCombat




def load_n_estimate(df, nnpc, mp, GRM, PC_effect="fixed", std=True, Method="AdjHE", random_groups=None, silent=False, homo=True, gcta=None, qcovar=None, covar_discrete=None, all_cols=None):
    """
    Estimates heritability, but solves a full OLS problem making it slower than the closed form solution. Takes 
    a dataframe, selects only the necessary columns (so that when we do complete cases it doesnt exclude too many samples)
    residualizes the phenotype, then documents the heritability, standard error and some computer usage metrics.

    Parameters
    ----------
    df : pandas dataframe
        dataframe containing phenotype, covariates, and principal components.
    nnpc : int
        number of pcs to include.
    mp : int
        which phenotype to estimate on.
    GRM : np array
        the GRM with missingness removed.
    PC_effect : str
        specifies whether PCs are treated as "fixed", "mixed", or "random".
    std : bool, optional
        specifying whether standardization happens before heritability estimation. The default is True.
    Method: str
        specify which method of estimation to use AdjHE, PredLMM, or MOM
        Default is AdjHE
    random_groups : string, optional
        Variable to control for as a random effect, if applicable
    qcovar : list, optional
        List of quantitative covariates
    covar_discrete : list, optional
        List of discrete covariates
    all_cols : list, optional
        All available columns from covariate file (for GCTA auto-detection)

    Returns
    -------
    pandas dataframe containing:
        - heritability estimate
        - standard error the estimate
        - the phenotype
        - the number of pcs included
        - The covariates included 
        - time for analysis
        - maximum memory usage
    """

    # Change nnpc to a number
    if nnpc == None :
        nnpc = 0
    
    # Get PC columns flexibly - match any case variation
    pc_cols = [c for c in df.columns if c.lower().startswith("pc")]
    # Remove duplicates while preserving order
    seen = set()
    pc_cols = [x for x in pc_cols if not (x in seen or seen.add(x))]
    # Sort numerically
    pc_cols = sort_pc_columns(pc_cols)
    logging.debug(f"PC columns found in df (deduplicated): {pc_cols[:20]} (total: {len(pc_cols)})")
    logging.debug(f"Requested nnpc: {nnpc}")
    if len(pc_cols) >= nnpc:
        pc_cols = pc_cols[:nnpc]
    else:
        # Use pc1, pc2 format to match actual column names in PC file
        pc_cols = ["pc" + str(p) for p in range(1, nnpc +1)]

    logging.debug(f"PC columns after slicing: {pc_cols}")

    # Remove duplicates while preserving order
    seen = set()
    pc_cols = [x for x in pc_cols if not (x in seen or seen.add(x))]

    # Build fixed_effects from qcovar and covar_discrete
    fixed_effects = []
    if qcovar:
        # Filter out PC columns (they're handled separately via nnpc)
        pc_in_qcovar = [c for c in qcovar if c.lower().startswith("pc")]
        non_pc_qcovar = [c for c in qcovar if not c.lower().startswith("pc")]
        if pc_in_qcovar:
            logging.warning(f"PC columns {pc_in_qcovar} in qcovar are ignored. Use 'npc' parameter to control PCs.")
        fixed_effects.extend(non_pc_qcovar)
    if covar_discrete:
        # Filter out PC columns
        pc_in_covar_discrete = [c for c in covar_discrete if c.lower().startswith("pc")]
        non_pc_covar_discrete = [c for c in covar_discrete if not c.lower().startswith("pc")]
        if pc_in_covar_discrete:
            logging.warning(f"PC columns {pc_in_covar_discrete} in covar_discrete are ignored. Use 'npc' parameter to control PCs.")
        fixed_effects.extend(non_pc_covar_discrete)

    # Remove duplicates while preserving order
    seen = set()
    fixed_effects = [x for x in fixed_effects if not (x in seen or seen.add(x))]

    # GCTA is the only method that doesn't need to specify whether npcs are fixed or mixed.
    if Method != "GCTA" :
        if PC_effect == "fixed" :
            # Add PC columns, avoiding duplicates
            for pc in pc_cols:
                if pc not in fixed_effects:
                    fixed_effects.append(pc)
            nnpc = 0
        if PC_effect == "mixed" :
            for pc in pc_cols:
                if pc not in fixed_effects:
                    fixed_effects.append(pc)
        # Random is implicitly handled

    # Create formula string
    if len(fixed_effects) != 0:
        RHS = " + ".join(fixed_effects)
    else :
        RHS = " 1"
    
    # Make formula
    form = f'{mp} ~ {RHS}'
    logging.debug("First Momment formula is " + form)
    
    # Select method of estimation
    try :
        if Method == "AdjHE":
            # AdjHE projects away covariates to start
            resid = smf.ols(formula=form, data=df, missing='drop').fit().resid
            resid.name = "resid"
            temp = df.merge(resid, left_index = True, right_index =True, how = "inner")
            temp[mp] = temp["resid"]
            nonmissing = df[df.IID.isin(temp.IID)].index
            GRM_nonmissing = GRM[nonmissing, :][:, nonmissing]
            result = AdjHE(A = GRM_nonmissing, df=temp, mp = mp, random_groups = random_groups, npc= nnpc, std=std)
            result["N"] = len(temp)

        elif (Method == "GCTA") or (Method == "HEreg"):
            from Estimate.estimators.GCTA_wrapper import gcta as gcta_path
            result = GCTA(df, nnpc, mp, GRM, gcta=gcta_path, method=Method, silent=False, qcovar=qcovar, covar_discrete=covar_discrete, all_available_cols=all_cols)
            result["N"] = len(df)
            
        elif Method == "SWD":
            # SWD projects away sites then projects away covaraites
            resid = smf.ols(formula= f'{mp} ~ {random_groups}', data=df, missing='drop').fit().resid
            resid.name = "resid"
            temp = df.merge(resid, left_index = True, right_index =True, how = "inner")
            temp[mp] = temp["resid"]
            
            resid = smf.ols(formula=form, data=temp, missing='drop').fit().resid 
            resid.name= "resid2"
            temp = temp.merge(resid, left_index= True, right_index = True,how = "inner")
            temp[mp] = temp["resid2"]
            nonmissing = df[df.IID.isin(temp.IID)].index
            GRM_nonmissing = GRM[nonmissing, :][:, nonmissing]

            result = AdjHE(A = GRM_nonmissing, df = temp, mp = mp, random_groups = None, npc=nnpc, std=False)
            result["N"] = len(temp)

        elif Method in ["Combat", "Covbat"]:
            # AdjHE projects away covariates to start
            result = AdjHE(A = GRM, df=df, mp = mp, random_groups = None, npc= nnpc, std=std)
            result["N"] = len(df)


        else:
            logging.error("Not an accepted method of estimation: " + Method)
            result = {}

        result["pheno"] = mp
        result["N"] = result.get("N", len(df))
        return pd.DataFrame(result, index=[0])

    except np.linalg.LinAlgError :
        logging.error("Singular Matrix")
        pass
    except TypeError :
        logging.error("Muffed estimate")
        pass

class h2Estimation():
    def __init__(self, args= None,  k=0, ids=None):
        if args == None:
            logging.info("Enter preloaded values...")
            self.df = None
            self.GRM = None
            self.phenotypes = "Y"
            self.simulation = True
            self.args = None

        else:
            logging.info("Loading data...")
            self.args = args
            self.df, self.GRM, self.phenotypes, self.ids = load_everything(
                args= self.args, k= k)
            self.simulation = False
            # Store all available columns from covariate file for GCTA (for auto-detection of discrete vs quantitative)
            self.all_covar_cols = None
            if args.get("covar") and args["covar"] != "None":
                covar_input = args["covar"]
                if isinstance(covar_input, str):
                    covar_files = [f.strip() for f in covar_input.split(',')]
                elif isinstance(covar_input, list):
                    covar_files = covar_input
                else:
                    covar_files = [covar_input]
                all_covar_cols = set()
                for cov_file in covar_files:
                    try:
                        import os
                        ext = os.path.splitext(cov_file)[1].lower()
                        if ext in ['.tsv', '.tab']:
                            sep = '\t'
                        elif ext == '.csv':
                            sep = ','
                        elif ext in ['.phe', '.pheno', '.phen', '.txt', '.covar']:
                            sep = None  # whitespace
                        else:
                            sep = None
                        cov_df = pd.read_csv(cov_file, sep=sep, nrows=1)
                        all_covar_cols.update(c for c in cov_df.columns if c not in ["FID", "IID"])
                    except:
                        pass
                self.all_covar_cols = list(all_covar_cols)
            #self.df = pd.merge(self.ids, self.df, on = ["FID", "IID"], how = "left")

    def estimate(self, PC_effect = "fixed"):
        args = self.args

        # Build fixed_effects from qcovar and covar_discrete
        qcovar = args.get("qcovar") or []
        covar_discrete = args.get("covar_discrete") or []
        fixed_effects = qcovar + covar_discrete
        # Remove duplicates while preserving order
        seen = set()
        fixed_effects = [x for x in fixed_effects if not (x in seen or seen.add(x))]

        # Create list of covariate sets to regress over
        if not fixed_effects:
            fixed_combos = [[]]
        else:
            # Create the sets of covariates over which we can loop
            # This will return a list of lists of covariate names to regress on
            fixed_combos = [fixed_effects[0:idx+1] for idx, c in enumerate(fixed_effects)]
            # If we don't want to loop, just grab the last item of the generated list assuming the user wants all of those variables included
            if not args["loop_covars"]:
                fixed_combos = [fixed_combos[-1]]

        if args["mpheno"] == "all":
            self.mpheno = self.phenotypes
        else:
            # make them lowercase
            self.mpheno = args["mpheno"] 
            # Convert numeric mpheno to column names if needed
            if self.mpheno and isinstance(self.mpheno[0], int):
                self.mpheno = ["pheno_" + str(m) for m in self.mpheno]
            
        logging.info("Estimating with " + args["Method"])
        
        if args["random_groups"] != "None":
            logging.info("RV: " + str(args["random_groups"]))
        # create empty list to store heritability estimates
        results = pd.DataFrame()

        # project GRM onto pc space
        # pcs = pd.DataFrame(PCA(n_components=20).fit_transform(self.GRM))
        # pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
        # # add the pcs to the dataframe
        # self.df = pd.concat([self.df, pcs], axis=1)

        # Forcing type to be integer for a little easier use
        if args["npc"] == None:
            npc = [0]
            
        # Adjust data if any Combat based method is wanted
        if args["Method"] in ["Combat", "Covbat"] :

            logging.info("Method: " + args["Method"])
            # Get PC columns flexibly
            all_pc_cols = [c for c in self.df.columns if c.lower().startswith("pc") and c.lower() not in ["fid", "iid"]]
            qcovar = args.get("qcovar") or []
            covar_discrete = args.get("covar_discrete") or []
            all_covars = qcovar + covar_discrete
            if all_pc_cols:
                FEs = all_pc_cols + all_covars if all_covars else all_pc_cols
            elif args["npc"]:
                FEs = ["pc" + str(i + 1) for i in range(max(args["npc"]))] + all_covars if all_covars else ["pc" + str(i + 1) for i in range(max(args["npc"]))]
            else:
                FEs = all_covars
            random_groups = args["random_groups"]
            no_missing = self.df[["FID", "IID"] + self.mpheno + FEs + [random_groups]].dropna()
            
            # Run neuroCombat: input shape (n_phenotypes, n_subjects), output same shape
            # After .T: shape is (n_subjects, n_phenotypes)
            transformed_data = neuroCombat(dat=no_missing[self.mpheno].T,
                                           covars=no_missing[FEs + [args["random_groups"]]],
                                           batch_col=args["random_groups"])["data"].T
            
            args["random_groups"] = None
            
            # Filter self.df to only rows that were used (no missing)
            nonmissing_idx = no_missing.index
            self.GRM = self.GRM[nonmissing_idx, :][:, nonmissing_idx]
            self.df = self.df.loc[nonmissing_idx].copy()
            
            # Assign transformed phenotype data back to dataframe
            # Make sure shapes match: transformed_data should be (n_subjects, n_phenotypes)
            if len(transformed_data.shape) == 1:
                transformed_data = transformed_data.reshape(-1, 1)
            
            # Ensure self.df has the phenotype columns before assignment
            # If transformed_data has different shape than expected, adjust
            n_subjects, n_phenotypes = transformed_data.shape
            if n_phenotypes == len(self.mpheno):
                self.df[self.mpheno] = transformed_data
            else:
                # Handle case where number of phenotypes changed (e.g., PCA was applied)
                logging.warning(f"Transformed data has {n_phenotypes} phenotypes, expected {len(self.mpheno)}")
                # For now, just use first n_phenotypes columns
                for i, pheno in enumerate(self.mpheno[:n_phenotypes]):
                    self.df[pheno] = transformed_data[:, i]
            
            # After transforming, Combat and Covbat procedures proceed just as the basic AdjHE estimator

        
        logging.info("Beginning estimation")
        # Calculate total iterations for progress bar
        total_iters = len(fixed_combos) * len(self.mpheno) * len(args["npc"])
        pbar = tqdm(total=total_iters, desc="Estimating", unit="phenotype")
        
        # Loop over each set of covariate combos
        for covs in fixed_combos:
            # For each set of covariates recalculate the projection matrix
            logging.debug(covs)
            # loop over all combinations of pcs and phenotypes
            for mp, nnpc in itertools.product(self.mpheno, args["npc"]):
                pbar.set_description(f"Estimating: {mp}")
                pbar.update(1)
                
                start_est = timeit.default_timer()


                
                try: 
                    C = "+".join(covs)
                except TypeError :
                    C = "None"
                    
                r = {"Pheno": mp, 
                          "PCs" : nnpc,
                          "Covariates" : C}

                if (not args["Naive"]) or (random_groups == None):
                    r = load_n_estimate(df=self.df, nnpc=nnpc,
                                        mp=mp, GRM=self.GRM, std=True, Method=args["Method"],
                                        random_groups=args["random_groups"], homo=True, PC_effect = PC_effect,
                                        qcovar=args.get("qcovar"), covar_discrete=args.get("covar_discrete"),
                                        all_cols=self.all_covar_cols)

                else:
                    # Empty results list
                    sub_results = pd.DataFrame({"h2": [],
                                                "Size": []})

                    # loop over  all sites
                    groups= np.unique(self.df[args["random_groups"]])
                    groups = groups[~np.isnan(groups)]
                    for group in tqdm(groups, desc=f"Analyzing {mp} across sites", leave=False):
                        try : 
                            # Grab the portion that lies within a given site
                            sub_df = self.df.loc[self.df[args["random_groups"]] == group, :].reset_index(drop = True)
                            # Get size
                            sub_n = sub_df.shape[0]
                            sub_GRM = self.GRM[self.df[args["random_groups"]] == group,:][:,self.df[args["random_groups"]] == group]
                            # Find PC's individually for each site
                            if nnpc != 0:
                                pcs = pd.DataFrame(PCA(n_components=30).fit_transform(np.asarray(sub_GRM)))
                                pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
                                # Drop previous pcs
                                sub_df = sub_df.loc[:,~sub_df.columns.str.startswith('pc_')]
                                #Add site specific PC's
                                sub_df = pd.concat([sub_df, pcs], axis=1)



                            # Estimate just on the supsample
                            sub_result = load_n_estimate(df=sub_df, nnpc=nnpc, mp=mp, GRM=sub_GRM, std=True, Method=args["Method"], random_groups=None,
                                                     silent=True, homo=True, PC_effect = PC_effect,
                                                     qcovar=args.get("qcovar"), covar_discrete=args.get("covar_discrete"),
                                                     all_cols=self.all_covar_cols)
                            sub_result = pd.DataFrame({"h2": [sub_result["h2"][0]],
                                                    "Size": [sub_n],
                                                    "N": [sub_result.get("N", sub_n)]})
                            # Add to the list of estimates
                            sub_results = pd.concat([sub_results, sub_result], axis=0)
                        except ValueError :
                            logging.error("Not estimated on this subgroups since there wasn't enough samples")

                    # Pool the estimates
                    sub_results["nh2"] = (
                        sub_results["Size"] * sub_results["h2"]) / self.GRM.shape[0]
                    h2 = np.sum(sub_results["nh2"])
                    r["h2"] = h2
                    r["var(h2)"] =  np.var(sub_results["h2"])
                    r = pd.DataFrame(r, index=[0])

                # Get memory for each step (in Mb) (This is a little sketchy)
                logging.debug("Started" + str(start_est)) 
                end_est = timeit.default_timer()
                logging.debug("Ended" + str(end_est))
                logging.debug(end_est - start_est)
                
                time = [end_est - start_est]
                mem = [0]  # Placeholder for memory usage (removed resource dependency)
                r["PCs"] = nnpc
                r["time"] = time
                r["mem"] = mem
                r["N"] = r.get("N", np.nan)
                results = pd.concat([results, r], ignore_index=True)

        pbar.close()
        self.results = results
        return results # , sub_results




