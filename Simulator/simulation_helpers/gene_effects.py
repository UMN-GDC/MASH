#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:42:57 2023

@author: christian
"""
import numpy as np
import matplotlib.pyplot as plt

def sim_gen_effects(rng, genotypes, alpha = -1):
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
    return Gene_contrib, causals, causal_eff







#%% different function

    # plot distribution of causal effects,
    # causal effect against index, and 
    # casual effect vs allele frequency 
    # in a three panel plot
    if vizualize:
        fig, axs = plt.subplots(2, 2)
        axs[0,0].hist(causal_eff, bins = 100)
        axs[0,0].set_title('Distribution of causal effects')
        axs[0,1].scatter(causals, causal_eff)
        axs[0,1].set_title('Causal effect vs index')
        axs[1,1].scatter(freqs, causal_eff)
        axs[1,1].set_title('Causal effect vs allele frequency')
    return None





