#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:42:57 2023

@author: christian
"""
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf


# Need to add shared here then select causal from shared and nonshared 
def sim_pheno(rng, genotypes, df, h2Hom, h2Het, alpha = -1, phenoname = "Y0", causals= None, prop_causal= [0.1, 0.1], sharedIdx = None):
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
    prop_causal : list, optional
        if causals are unspecified, then specify the number of SNPs to select as causal from shared and unshared regions. The default is 0.1.
    clusters_differ : bool, optional
        do clusters have different heritabilities? The default is False.
    maf_filter : float, optional
        filter out SNPs with MAF < maf_filter and MAF > 1-maf_filter. This is only used when user did not specify causals explicitly. The default is 0.1.
    vizualize : boolean, optional
        plot the distribution of causal effects. The default is False.
    sharedIdx: boolean array 
        Index of SNPs that have shared frequencies between populations
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
        nSNPs = genotypes.shape[1]
        causals = np.repeat(False, nSNPs)
        causals[:int(prop_causal[0] * sum(sharedIdx))] = True
        causals[int(sum(sharedIdx)): int(sum(sharedIdx) + ((nSNPs - sum(sharedIdx)) * prop_causal[1]))] = True

        nCausal = sum(causals)
    
    Xcausal = np.matrix(genotypes[:, causals])
    freqs = np.asarray(np.mean(Xcausal, axis = 0)/2).flatten()

    # sample shared (homogeneous) effects from normal with variance proportional to some function of global allele frequency
    prop = (freqs * (1-freqs))**alpha
    homo_eff = rng.normal(np.repeat(0, nCausal), np.sqrt(h2Hom / nCausal * (2 * prop) ** alpha  * 2 ** (alpha+ 1)), size =  nCausal)
    # make sure no infinities
    homo_eff[np.isinf(homo_eff)] = 0
    homo_contrib = np.array(np.dot(Xcausal, homo_eff)).flatten()
    df[f"homo_contrib{phenoname}"] = homo_contrib * np.sqrt(h2Hom / np.var(homo_contrib))
    nclusts = df.subj_ancestries.nunique()
    cluster_eff = np.zeros((nCausal, nclusts))
    cluster_contrib = np.zeros(nsubjects)
    errors = np.zeros(nsubjects)
    # cluster specific effects
    if nclusts > 1 :
        for cluster in range(nclusts) :
            cluster_position = df.subj_ancestries == cluster 
            cluster_freq = np.asarray(np.mean(genotypes[cluster_position,:][:,causals], axis = 0)/2).flatten()
            prop = (cluster_freq * (1- cluster_freq)) ** alpha
            cluster_eff[:,cluster] = rng.normal(np.repeat(0, nCausal), prop, size =  nCausal)
            cluster_eff[np.isinf(cluster_eff)] = 0
            cluster_contrib[cluster_position] = np.array(np.dot(Xcausal[cluster_position, :], cluster_eff[:,cluster])).flatten()
            # rescale
            cluster_contrib[cluster_position]= cluster_contrib[cluster_position] * np.sqrt(h2Het[cluster]/ np.var(cluster_contrib[cluster_position]))

            cluster_error = rng.normal(0, 1, cluster_position.sum())
            errors[cluster_position] = cluster_error * np.sqrt((1- h2Hom - h2Het[cluster]) / np.var(cluster_error))
    else :
        errors = rng.normal(0, np.sqrt(1-h2Hom), nsubjects)
        # multiply errors such that its variance is equal to 1-h2Hom
    if h2Het[0] == 0  :
        cluster_contrib =0
    df[f"cluster_contrib{phenoname}"] = cluster_contrib
    df[f"errors{phenoname}"] = errors
    df["Xc"] = rng.uniform(0, 1, nsubjects)
    Beta_c = 0.5
    df[f"Covar_contrib{phenoname}"] = Beta_c * df["Xc"]
    
    # sum across all columns with phenoname at the end of the column name using regex
    df[str(phenoname)] = df.filter(regex = f"{phenoname}").sum(axis = 1)
    
    df[str(phenoname)] -= df[str(phenoname)].mean()
    

    return df, causals
