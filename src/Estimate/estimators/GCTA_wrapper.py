#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 07:35:15 2022

@author: christian
"""
import sys
import os
import logging
import subprocess
import numpy as np
import pandas as pd 
import re
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")


def sort_pc_columns(pc_cols):
    """Sort PC column names numerically (pc1, pc2, ..., pc10, pc11)."""
    def get_pc_num(col):
        match = re.search(r'pc(\d+)', col.lower())
        return int(match.group(1)) if match else float('inf')
    return sorted(pc_cols, key=get_pc_num)


# Find GCTA
mash_gcta_paths = [
    "/projects/standard/gdc/public/envs/MASH/bin/gcta64",
    "/projects/standard/gdc/public/envs/MASH/bin/gcta",
    os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin", "gcta64"),
    os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin", "gcta"),
]
for path in mash_gcta_paths:
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
        gcta = path
        break
else:
    try:
        gcta = subprocess.run(["which", "gcta64"], capture_output=True
                           ).stdout.split()[1].decode("utf-8") 
    except IndexError : 
        try :
            gcta = subprocess.run(["which", "gcta"], capture_output=True
                                  ).stdout.split()[1].decode("utf-8") 
        except IndexError:
            try :
                gcta = subprocess.run(["which", "gcta"], capture_output=True
                                      ).stdout.decode("utf-8") 
            except IndexError:
                logging.error("GCTA was not found. Please install GCTA and make sure it is in your path")
                sys.exit(1)

# Expand tilde to full path
gcta = os.path.expanduser(gcta)




#%%
def GCTA(df, nnpc, mp, GRM, gcta, method = "GCTA", silent=False, qcovar=None, covar_discrete=None, all_available_cols=None):
    import tempfile
    import hashlib

    # Use all available columns from covariate file for auto-detection if not provided
    if all_available_cols is None:
        all_available_cols = list(set(df.columns) - {"FID", "IID"})
    available_cols = set(all_available_cols)

    # Validate that specified qcovar and covar_discrete columns exist
    if qcovar is not None:
        qcovar = [c for c in qcovar if c in available_cols]
        missing = set(qcovar) - available_cols
        if missing:
            logging.warning(f"qcovar columns not found in data: {missing}")
    if covar_discrete is not None:
        covar_discrete = [c for c in covar_discrete if c in available_cols]
        missing = set(covar_discrete) - available_cols
        if missing:
            logging.warning(f"covar_discrete columns not found in data: {missing}")

    # Create a unique temp directory based on hash of inputs to avoid collisions
    input_hash = hashlib.md5(f"{mp}{nnpc}{qcovar}{covar_discrete}".encode()).hexdigest()[:8]
    temp_dir = tempfile.mkdtemp(prefix=f"gcta_{input_hash}_")
    temp_name = os.path.join(temp_dir, "gcta_temp")

    # write the phenotype file
    df[["FID", "IID", mp]].to_csv(temp_name + "_pheno.txt", sep = " ", header = False, index= False, na_rep = "NA")
    
    # Store copy of df for covariate file writing (before truncation)
    df_for_covars = df.copy()
    
    # Pre-compute discrete classification BEFORE truncating df
    # This is needed because we need access to all covariate columns for classification
    col_uniqueness = {}
    all_covar_cols_in_df = [c for c in all_available_cols if c in df.columns]
    for col in all_covar_cols_in_df:
        col_uniqueness[col] = len(df[col].unique())
    
    # Select the remaining variables of interest
    # Match PCs flexibly - look for columns starting with "pc" or "PC" (case-insensitive)
    pc_cols = [c for c in df.columns if c.lower().startswith("pc")]
    pc_cols = sort_pc_columns(pc_cols)
    if nnpc > 0 and len(pc_cols) >= nnpc:
        pcs = pc_cols[:nnpc]
    else:
        # Use pc1, pc2, etc. format to match the actual column names
        pcs = ["pc" + str(s + 1) for s in range(nnpc)]

    # Build list of all covariate columns to keep
    all_covar_cols = []
    if qcovar:
        all_covar_cols.extend(qcovar)
    if covar_discrete:
        all_covar_cols.extend(covar_discrete)

    # Store all covariate columns before truncating df (for writing covariate files later)
    all_covars_in_df = all_covar_cols_in_df
    df = df[["FID", "IID"] + all_covar_cols + pcs]

    # Determine qcovar and covar_discrete columns: explicit > inferred
    # Only auto-detect from columns NOT already included in all_covar_cols + pcs (avoid duplicates)
    already_used = set(all_covar_cols) | set(pcs)
    remaining_cols = [c for c in all_available_cols if c not in already_used]

    if qcovar is not None or covar_discrete is not None:
        if qcovar is None:
            qcovar = []
        if covar_discrete is None:
            covar_discrete = []
    else:
        discrete = []
        for col in remaining_cols:
            n_unique = col_uniqueness.get(col, 0)
            discrete.append((n_unique < 35) and (n_unique > 1))

        cont = [not v for v in discrete]

        qcovar_cols = [col for i, col in enumerate(remaining_cols) if cont[i] and col not in ["FID", "IID"]]
        covar_cols = [col for i, col in enumerate(remaining_cols) if discrete[i] and col not in ["FID", "IID"]]
        qcovar = qcovar_cols
        covar_discrete = covar_cols

    # Write covariate files if there are any in either category
    # Use df_for_covars which contains all covariate columns
    covar_discrete_in_df = [c for c in covar_discrete if c in all_covars_in_df]
    qcovar_in_df = [c for c in qcovar if c in all_covars_in_df]

    # Include PCs as quantitative covariates for GCTA
    if nnpc > 0:
        for pc in pcs:
            if pc not in qcovar_in_df:
                qcovar_in_df.append(pc)

    if covar_discrete_in_df:
        df_for_covars[["FID", "IID"] + covar_discrete_in_df].to_csv(temp_name + "_Discrete.txt", sep = " ", header= False, index= False, na_rep = "NA")
    if qcovar_in_df:
        df_for_covars[["FID", "IID"] + qcovar_in_df].to_csv(temp_name + "_Cont.txt", sep = " ", header= False, index= False, na_rep = "NA")
    
    #######################
    # Write GRM and ids
    # Write IDs
    df[["FID", "IID"]].to_csv(temp_name + ".grm.id", sep = " ", header= False, index= False)
    
    # Write GRM to binary 
    (GRM[np.tril_indices(df.shape[0])]
     .astype("f4") # Relatedness is stored as a float of size 4 in the binary file 
     .tofile(temp_name + ".grm.bin")
    )
    ##############################
        
    # Format string for controlling variables
    covs = " "
    if qcovar_in_df:
        covs += " --qcovar " + temp_name + "_Cont.txt "
    if covar_discrete_in_df:
        covs += " --covar " + temp_name + "_Discrete.txt "
    
    if method == "GCTA":
        estimator = " --reml --reml-priors 0.05 0.95"
    else :
        estimator = " --HEreg" 
    
    # Combine all covariates for the result string
    all_covariates = []
    if qcovar_in_df:
        all_covariates.extend(qcovar_in_df)
    if covar_discrete_in_df:
        all_covariates.extend(covar_discrete_in_df)

    # run gcta
    bashcommand = f"{gcta} --grm {temp_name} --pheno {temp_name}_pheno.txt --mpheno 1 {estimator} --out {temp_name} {covs}"
    print(bashcommand)
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    __output, __error = process.communicate()

    # parse output for estimate
    try :
        if method == "GCTA" :
            df = pd.read_table(temp_name + ".hsq", sep="\t")
            result = {"h2": df.Variance[3], "var(h2)": df.SE[3]**2, "G" : df.Variance[0], "E" : df.Variance[1], "pval" :df.SE[8], "pheno": mp, "PCs": nnpc, "Covariates": "+".join(all_covariates), "time": np.nan,
                      "mem": np.nan, "N" : df.Variance[9]}
        if method == "HEreg" :
            df = pd.read_table(temp_name + ".HEreg", nrows=2, skiprows=1, sep="\s+")
            result = {"h2": df.Estimate[1], "var(h2)": df.SE_OLS[1]**2, "G" : np.nan, "E" : np.nan, "pval" :df.P_OLS[1], "pheno": mp, "PCs": nnpc, "Covariates": "+".join(all_covariates), "time": np.nan,
                 "mem": np.nan, "N": np.nan}

    except FileNotFoundError:
        try :
            if method == "GCTA" :
                logging.error("Estimations were not made. Trying again unconstrained and with seeded reml estimates")
                bashcommand = f"{gcta} --grm {temp_name} --pheno {temp_name}_pheno.txt --mpheno 1 {estimator} --out {temp_name} {covs} --reml-no-constrain --reml-maxit 200 --reml-priors 0.025 0.975"
                print(bashcommand)
                process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
                __output, __error = process.communicate()
                df = pd.read_table(temp_name + ".hsq", sep="\t")
                result = {"h2": df.Variance[3], "var(h2)": df.SE[3]**2, "G" : df.Variance[0], "E" : df.Variance[1], "pval" :df.SE[8], "pheno": mp, "PCs": nnpc, "Covariates": "+".join(all_covariates), "time": np.nan,
                     "mem": np.nan, "N": df.Variance[9]}
        
        except FileNotFoundError:
            logging.error("Estimations were not made. Usually this is due to small sample sizes for GCTA")
            result = {"h2": np.nan, "var(h2)": np.nan, "pheno": mp, "PCs": nnpc, "Covariates": "+".join(all_covariates), "time": np.nan,
             "mem": np.nan, "N": np.nan}

    
    
    
    # tidy up by removing temporary directory
    import shutil
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    # Return the fit results
    return pd.DataFrame(result, index = [0])

        



