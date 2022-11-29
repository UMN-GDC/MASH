#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 07:35:15 2022

@author: christian
"""
import os
import subprocess
import pandas as pd 
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")

#%%

def GCTA(grm, pheno_file, cov_file, PC_file, covars, nnpc, mp) : 
    """
    A wrapper for the GCTA GREML heritability estimator.

    Parameters
    ----------
    grm : string
        file path prefix to GRM files.
    pheno_file : string
        file path to phenotype file.
    cov_file : string
        file path to covariates file.
    PC_file : string
        file path to covariates file.
    covars : list of strings
        list of covaraites to control for as fixed effects.
    nnpc : int
        number of pc's to control for as fixed effects.
    mp : string
        name of phenotype to estimate on.

    Returns
    -------
    pandas dataframe.
        Dataframe of estimated heritability and standard error

    """
    # Load and select covariates
    eigs = pd.read_csv(PC_file, sep = " ", header = None)
    cols = eigs.columns.tolist()[2:]
    eigs.columns = ["fid", "iid"] + [i -1  for i in cols]
    eigs = eigs.iloc[:,:nnpc+2]
    cov = pd.read_csv(cov_file, sep = " ", header= 0)
    cov.columns = [col.lower() for col in cov.columns.tolist()]
    cov = cov[["fid", "iid"] + covars]
    cov = cov.merge(eigs, on = ["fid", "iid"])
    # Decide which are qcovars and which are numeric covars
    discrete = [len(cov[col].unique()) < 50 for col in cov]
    # Include FID IID
    discrete[0:2] = [True, True]
    cont = [not v for v in discrete]
    cont[0:2] = [True, True]
    
    # Svae temp files
    cov.iloc[:,discrete].to_csv("temp_Discrete.txt", sep = " ", header = False, index= False)
    cov.iloc[:,cont].to_csv("temp_Cont.txt", sep = " ", header= False, index= False)
    # Find GCTA
    gcta = "whereis gcta64"
    gcta, __ = subprocess.Popen(
        gcta.split(), stdout=subprocess.PIPE).communicate()
    gcta= gcta.split()[1].decode("utf-8")
    # run gcta
    out = "test"
    bashcommand = gcta + f" --grm {grm} --pheno {pheno_file} --mpheno 1 --reml --out {out} --qcovar temp_Cont.txt --covar temp_Discrete.txt"
    process = subprocess.Popen(bashcommand.split(), stdout=subprocess.PIPE)
    __output, __error = process.communicate()

    #
    df = pd.read_table(out + ".hsq", sep="\t").query( "Source == 'V(G)/Vp'").reset_index()
    
    result = {"h2" : df.Variance[0],
              "SE" : df.SE[0],
              "Pheno" : mp,
              "PCs" : nnpc,
              "Covariates" : "+".join(covars),
              "Time for analysis(s)" : 0,
              "Memory Usage" : 0}
    
    # tidy up by removing temporary files
    os.remove("temp_Discrete.txt")
    os.remove("temp_Cont.txt")
    
    # Return the fit results
    return(pd.DataFrame(result, index = [0]))




