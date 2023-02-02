#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 07:35:15 2022

@author: christian
"""
import sys
import os
import subprocess
import numpy as np
import pandas as pd 
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")


# Find GCTA
gcta = subprocess.run(["whereis", "gcta64"], capture_output=True
                      ).stdout.split(
                          )[1].decode("utf-8")


#%%
def GCTA(df, covars, nnpc, mp, GRM, gcta, silent=False):
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
    # Specify information about binary GRM format
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
    
    
    # Write IDs
    df[["FID", "IID"]].to_csv(temp_name + ".grm.id", sep = " ", header= False, index= False)
    n = df.shape[0]
    
    l= np.tril_indices(n)
    
    # Write GRM to binary 
    GRM[l].astype("f4").tofile(temp_name + ".grm.bin")    
    ##############################
        
    # Format string for controlling variables
    covars = " "
    if os.path.exists(temp_name + "_Cont.txt") :
        covars += " --qcovar " + temp_name + "_Cont.txt "
    if os.path.exists(temp_name + "_Discrete.txt") : 
        covars += " --covar " + temp_name + "_Discrete.txt "
    
    
    # run gcta
    bashcommand = gcta + " --grm " + temp_name + " --pheno " + temp_name + "_pheno.txt --mpheno 1 --reml --out " + temp_name + " " + covars
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    __output, __error = process.communicate()

    # parse output for estimate
    try :
        df = pd.read_table(temp_name + ".hsq", sep="\t").query( "Source == 'V(G)/Vp'").reset_index()
        result = {"h2" : df.Variance[0],
                  "var(h2)" : df.SE[0]**2,
                  "ss" : 0}

    except FileNotFoundError:
        print("Estimations were not made. Usually this is due to small sample sizes for GCTA")
        result = {"h2": np.nan, "SE": np.nan, "Pheno": mp, "PCs": nnpc, "Covariates": "+".join(covars), "Time for analysis(s)": np.nan,
             "Memory Usage": np.nan}

    
    
    
    # tidy up by removing temporary files
    if os.path.exists(temp_name + "_Discrete.txt") : 
        os.remove(temp_name + "_Discrete.txt")
    if os.path.exists(temp_name + "_Cont.txt") : 
        os.remove(temp_name + "_Cont.txt")
    if os.path.exists(temp_name + ".hsq") : 
        os.remove(temp_name + ".hsq")
    if os.path.exists(temp_name + ".log") : 
        os.remove(temp_name + ".log")
    if os.path.exists(temp_name + "_pheno.txt") :
        os.remove(temp_name + "_pheno.txt")
    if os.path.exists(temp_name + ".grm.bin") :
        os.remove(temp_name + ".grm.bin")
    if os.path.exists(temp_name + ".grm.id") :
        os.remove(temp_name + ".grm.id")


    
    # Return the fit results
    return pd.DataFrame(result, index = [0])

        



