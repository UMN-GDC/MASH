#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:42:57 2023

@author: christian
"""
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

def SNPEffectSampler(rng, genotypes, variance=0.5) :
    """
    Given a matrix of CAUSAL SNPs, generate a set of SNP effects that yields the desired variance (h2)
    for the contribution of G \cdot beta. 

    Parameters
    ----------
    rng : random number generator
        numpy random number generator.
    genotypes : array
        ( x ) array of unstandardized causal SNPs.
    variance : float, optional
        desired heritability. The default is 0.5.
    """
    freqs = genotypes.mean(axis = 0)/2
    beta = rng.normal(0, np.sqrt(variance/genotypes.shape[1] * (2*freqs* (1-freqs)) **(-1)))
    beta[np.isinf(beta)] = 0
    effect = np.dot(genotypes, beta)
    sig = np.var(effect)
    errors = rng.normal(0, np.sqrt(sig * (1- variance) / variance), effect.shape)
    # rescale to be precisely the variance we want
    effect = effect / np.sqrt(sig) * np.sqrt(variance)
    errors = errors / errors.std() * np.sqrt(1-variance)
    return beta, effect, errors


def sim_pheno(rng, genotypes, df, h2Hom, h2Het, phenoname = "Y0", causals= None):
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
   
    # sample causal SNPs if user didn't specify their own specific string of causal SNPs
    if causals is None :
        prop_causal = 0.1
        nCausal = int(nSNPs * prop_causal)
        causals = rng.choice(nSNPs, nCausal, replace=False, shuffle=False)
    else : 
        nCausal = len(causals)
    homo_eff, df["homo_contrib"], errors =  SNPEffectSampler(rng, genotypes[:,causals], variance=h2Hom) 

    nclusts = df.subj_ancestries.nunique()
    cluster_eff = np.zeros((nCausal, nclusts))
    cluster_contrib = np.zeros(nsubjects)
    # cluster specific effects
    if nclusts > 1 :
        for cluster in range(nclusts) :
            cluster_position = df.subj_ancestries == cluster 
            cluster_eff[:,cluster], cluster_contrib[cluster_position], errors[cluster_position] = SNPEffectSampler(rng,
                                                                                                                   genotypes[cluster_position,:][:,causals],
                                                                                                                   variance = h2Het[cluster]) 

    df["cluster_contrib"] = cluster_contrib
    df["errors"] = errors
    df["Xc"] = rng.uniform(0, 1, nsubjects)
    Beta_c = 1.5
    df["Covar_contrib"] = Beta_c * df["Xc"]
    df[["Covar_contrib", "homo_contrib", "cluster_contrib", "errors"]] =df[["Covar_contrib", "homo_contrib", "cluster_contrib", "errors"]] 
    phen = df[["Covar_contrib", "homo_contrib", "cluster_contrib", "errors"]].sum(axis = 1)
    df[str(phenoname)] = phen - phen.mean() 
    return df, causals, homo_eff, cluster_eff 

