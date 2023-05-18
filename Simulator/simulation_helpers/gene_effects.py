#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:42:57 2023

@author: christian
"""
import numpy as np


def sim_gen_effects(rng, genotypes, causals = [], prop_causal=0.1, variance_propto_frequency = False, maf_filter = 0.1, clusters_differ = False):
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
    Returns
    -------
    Gene_contrib : array 
        1-D array containing the contribution of genetics to each subjects phenotype.
    causals : array
        1-D array of the genotype index that was selected to be causal
    """
    nsubjects, nSNPs = genotypes.shape

    # randomly sample causals only if 
    if causals == []:
        nCausal = int(nSNPs * prop_causal)
        # select causal snos
        causals = rng.choice(nSNPs, nCausal, replace=False, shuffle=False)
        Xcausal = np.matrix(genotypes[:, causals])
        # calculate frequencies
        freqs = np.asarray(np.mean(Xcausal, axis = 0)/2).flatten()
        # filter Xcausal by freqs that are greater than maf_filter, and less than 1-maf_filter
        Xcausal = Xcausal[:, (freqs > maf_filter) & (freqs < 1-maf_filter)]
        nCausal = Xcausal.shape[1]

    else :
        Xcausal = np.matrix(genotypes[:, causals])
        nCausal = Xcausal.shape[1]

    freqs = np.asarray(np.mean(Xcausal, axis = 0)/2)


    # sample effects from normal with variance proportional to some function of allele frequency
    if variance_propto_frequency : 
        alpha = -1
        prop = (freqs * (1-freqs))**alpha
        prop = prop.reshape(prop.shape[1])
        causal_eff = rng.normal(np.repeat(0, Xcausal.shape[1]), prop, size =  Xcausal.shape[1])
        # make sure no infinities
        causal_eff[np.isinf(causal_eff)] = 0
        Gene_contrib = np.array(np.dot(Xcausal, causal_eff)).flatten()


    else:
        # sim effect from each SNP
        causal_eff = rng.normal(0, 1, (Xcausal.shape[1], 1))
        Gene_contrib = np.array(Xcausal * causal_eff).flatten()

    return Gene_contrib, causals

