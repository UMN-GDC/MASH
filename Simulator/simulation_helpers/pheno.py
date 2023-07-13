#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:42:57 2023

@author: christian
"""
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

def sim_gen_effects(rng, genotypes, df, alpha = -1):
    """
    

    Parameters
    ----------
    rng : random number generatorng
        numpy random number generator.
    genotypes : array
        ( x ) array of unstandardized genotypes.
    causals : list, optional
        list of causal snp names.
    df : dataframe
        dataframe containing all covariates and phenotypes.
    alpha : float, optional
        exponent for the dependency between SNP frequency and effect size. The default is -1
    prop_causal : float, optional
        if causals are unspecified, then specify the number of SNPs to randomly select as causal. The default is 0.1.
    clusters_differ : bool, optional
        do clusters have different heritabilities? The default is False.
    maf_filter : float, optional
        filter out SNPs with MAF < maf_filter and MAF > 1-maf_filter. This is only used when user did not specify causals explicitly. The default is 0.1.
    vizualize : boolean, optional
        plot the distribution of causal effects. The default is False.
    Returns
    -------
    Gene_contrib : array 
        1-D array containing the contribution of genetics to each subjects phenotype.
    causals : array
        1-D array of the genotype index that was selected to be causal
    snp_effects: array
        1-D array of SNP effects
    """
    nsubjects, nSNPs = genotypes.shape
    prop_causal = 0.1
    # randomly sample causals only if 
    nCausal = int(nSNPs * prop_causal)
    # select causal snos
    causals = rng.choice(nSNPs, nCausal, replace=False, shuffle=False)
    Xcausal = np.matrix(genotypes[:, causals])
    # calculate frequencies
    freqs = np.asarray(np.mean(Xcausal, axis = 0)/2).flatten()

    # sample effects from normal with variance proportional to some function of allele frequency
    prop = (freqs * (1-freqs))**alpha
    causal_eff = rng.normal(np.repeat(0, Xcausal.shape[1]), prop, size =  Xcausal.shape[1])
    # make sure no infinities
    causal_eff[np.isinf(causal_eff)] = 0
    Gene_contrib = np.array(np.dot(Xcausal, causal_eff)).flatten()
    nclusts = df.subj_ancestries.nunique()
    print(nclusts)
    # Regress out the specified pcs from the Gene_contrib
    if nclusts > 1 :
        pcs = " + ".join(["pc_" + str(i) for i in range(1, int(np.floor((nclusts+1 ) ** 0.75)))])
        print(pcs)
        mod = smf.ols(formula = 'Gene_contrib ~ ' + pcs, data = df).fit()
        # Get the residuals
        resid_eff = mod.resid
        # Get the predicted values
        PC_eff = mod.predict()
    else : 
        resid_eff = Gene_contrib
        PC_eff = 0 
    return PC_eff, resid_eff, causals, causal_eff

def sim_pheno(rng, df,  h2, phenoname = "Y0")  : 
    rng = np.random.default_rng(rng)
    # Scale df resid_eff so that it has a variance of h2
    df["resid_eff"] = df['resid_eff'] * np.sqrt(h2/ np.var(df["resid_eff"]))
    df["errors"] = rng.normal(0, 1, df.shape[0])
    df["errors"] = df["errors"] * np.sqrt((1-h2) / np.var(df["errors"]))
    # df["PC_eff"] = df["PC_eff"] * np.sqrt(h2 / np.var(df["PC_eff"])) / 2
    df[str(phenoname)] = df[["Covar_contrib", "PC_eff", "resid_eff", "errors"]].sum(axis = 1)


    return df[["PC_eff", "resid_eff", "errors", str(phenoname)]] 
