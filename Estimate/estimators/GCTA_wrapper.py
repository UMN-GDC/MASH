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
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")


# Find GCTA
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




#%%
def GCTA(df, covars, nnpc, mp, GRM, gcta, method = "GCTA", silent=False):
    # Store random integer to save to temp file name
    rng = np.random.default_rng()    
    numb = rng.integers(10000)
    temp_name = "temp" + str(numb)

    # write the phenotype file
    df[["FID", "IID", mp]].to_csv(temp_name + "_pheno.txt", sep = " ", header = False, index= False, na_rep = "NA")
    
    # Select the remaining variables of interest
    pcs =  ["pc_" + str(s + 1) for s in range(nnpc)]
    df = df[["FID", "IID"] + covars + pcs]


    # Decide which are qcovars and which are discrete covars, also elimnate completely in common variables
    discrete = [(len(df[col].unique()) < 35) and (len(df[col].unique()) > 1) for col in df]
    # Include FID IID
    cont = [not v for v in discrete]
    discrete[0:2] = [True, True]

    
    # Svae temp files if there were any covariates in either category
    if sum(discrete) > 2:
        df.iloc[:,discrete].to_csv(temp_name + "_Discrete.txt", sep = " ", header = False, index= False, na_rep = "NA")
    if sum(cont) > 2:
        df.iloc[:,cont].to_csv(temp_name + "_Cont.txt", sep = " ", header= False, index= False, na_rep = "NA")
    
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
    if os.path.exists(temp_name + "_Cont.txt") :
        covs += " --qcovar " + temp_name + "_Cont.txt "
    if os.path.exists(temp_name + "_Discrete.txt") : 
        covs += " --covar " + temp_name + "_Discrete.txt "
    
    if method == "GCTA":
        estimator = " --reml --reml-maxit 500 --reml-priors 0.1 0.9"
    else :
        estimator = " --HEreg" 
    
    # run gcta
    gcta = "/users/4/coffm049/software/bin/gcta64"
    bashcommand = f"{gcta} --grm {temp_name} --pheno {temp_name}_pheno.txt --mpheno 1 {estimator} --out {temp_name} {covs}"
    print(bashcommand)
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    __output, __error = process.communicate()

    # parse output for estimate
    try :
        if method == "GCTA" :
            df = pd.read_table(temp_name + ".hsq", sep="\t")
            result = {"h2": df.Variance[3], "var(h2)": df.SE[3]**2, "G" : df.Variance[0], "E" : df.Variance[1], "pval" :df.SE[8], "pheno": mp, "PCs": nnpc, "Covariates": "+".join(covars), "time": np.nan,
                 "mem": np.nan}
        if method == "HEreg" :
            df = pd.read_table(temp_name + ".HEreg", nrows=2, skiprows=1, sep="\s+")
            result = {"h2": df.Estimate[1], "var(h2)": df.SE_OLS[1]**2, "G" : np.nan, "E" : np.nan, "pval" :df.P_OLS[1], "pheno": mp, "PCs": nnpc, "Covariates": "+".join(covars), "time": np.nan,
                 "mem": np.nan}

    except FileNotFoundError:
        try :
            if method == "GCTA" :
                logging.error("Estimations were not made. Trying again unconstrained and with seeded reml estimates")
                bashcommand = f"{gcta} --grm {temp_name} --pheno {temp_name}_pheno.txt --mpheno 1 {estimator} --out {temp_name} {covs} --reml-no-constrain --reml-maxit 200 --reml-priors 0.025 0.975"
                print(bashcommand)
                process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
                __output, __error = process.communicate()
                df = pd.read_table(temp_name + ".hsq", sep="\t")
                result = {"h2": df.Variance[3], "var(h2)": df.SE[3]**2, "G" : df.Variance[0], "E" : df.Variance[1], "pval" :df.SE[8], "pheno": mp, "PCs": nnpc, "Covariates": "+".join(covars), "time": np.nan,
                     "mem": np.nan}
        
        except FileNotFoundError:
            logging.error("Estimations were not made. Usually this is due to small sample sizes for GCTA")
            result = {"h2": np.nan, "var(h2)": np.nan, "pheno": mp, "PCs": nnpc, "Covariates": "+".join(covars), "time": np.nan,
             "mem": np.nan}

    
    
    
    # tidy up by removing temporary files
    if os.path.exists(temp_name + "_Discrete.txt") : 
        os.remove(temp_name + "_Discrete.txt")
    if os.path.exists(temp_name + "_Cont.txt") : 
        os.remove(temp_name + "_Cont.txt")
    if os.path.exists(temp_name + "_pheno.txt") :
        os.remove(temp_name + "_pheno.txt")
    if os.path.exists(temp_name + ".grm.bin") :
        os.remove(temp_name + ".grm.bin")
    if os.path.exists(temp_name + ".grm.id") :
        os.remove(temp_name + ".grm.id")
    if os.path.exists(temp_name + ".HEreg") : 
        os.remove(temp_name + ".HEreg")
    if os.path.exists(temp_name + ".hsq") : 
        os.remove(temp_name + ".hsq")
    if os.path.exists(temp_name + ".log") : 
        os.remove(temp_name + ".log")


    
    # Return the fit results
    return pd.DataFrame(result, index = [0])

        



