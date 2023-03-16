#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def sim_genos(rng, ancestral_frequencies, cluster_frequencies, subject_ancestries, clusters_differ = False, prop_causal=0.1, maf_filter = 0.05):
    """
    Simulate genotypes of the subjects given their cluster ID and the respective cluster allele frequencies

    Parameters
    ----------
    rng : numpy random number generator Generator object
        random number generator.
    cluster_frequencies : numpy array
        (nsubjects x nclusts) array contianing allele frequencies for each cluster.
    subject_ancestries : pandas series
        vector containng which ancestry each subject will be sampled from.
    clusters_differ : bool, optional
        should clusters have different means or should all the variations in the phenotype take place in non-cluster informative regions.
        The default is False.
    prop_causal : float, optional
        between 0 and 1 specifyin the number of SNPs that are causal. The default is 0.1.
    maf_filter: float, optional
        between 0 and 0.5 specifyin the minimum allowable minor allele frequency. The default is 0.05.


    Returns
    -------
    tuple
    genotypes : numpy array
        (nsubjects x nSNPs array) containing subject genotypes.
    GRM : numpy array
        (nsubjects x nsubjects) array containing standardized genetic relatedness matrix .
    pcs : pandas dataframe
        dataframe containing loadings along all eigenvectors.
    causal_snps : numpy array
        numpy array of the causal SNPs

    """
    nSNPs = cluster_frequencies.shape[1]
    nclusts = cluster_frequencies.shape[0]
    nsubjects = subject_ancestries.shape[0]
    
    if nclusts == 1:
        # simulate genotypes
        genotypes = rng.binomial(n=np.repeat(2, nSNPs), p=cluster_frequencies[0],
                                 size=(nsubjects, nSNPs))
    elif clusters_differ:
        # simulate genotypes
        genotypes = rng.binomial(n=np.repeat(2, nSNPs), p= cluster_frequencies[subject_ancestries],
                                 size=(nsubjects, nSNPs))
        
    elif not clusters_differ: 
        anc_region = range(int(nSNPs *(1-prop_causal)))
        causal_region = range(max(anc_region)+1, nSNPs)
        anc_geno = rng.binomial(n=np.repeat(2, len(anc_region)),
                                 p= cluster_frequencies[subject_ancestries][:,anc_region],
                                 size=(nsubjects, len(anc_region)))
        causal_geno = rng.binomial(n=np.repeat(2, len(causal_region)),
                                 p=ancestral_frequencies[causal_region],
                                 size=(nsubjects, len(causal_region)))
        genotypes= np.concatenate((anc_geno, causal_geno), axis = 1)
        
        
    # keep SNPs with MAF greater than 0.05
    maf_filter = np.logical_and((np.sum(genotypes, axis=0) / (2 * nsubjects)) > maf_filter,
                                (np.sum(genotypes, axis=0) / (2 * nsubjects)) < 1-maf_filter)
    
    if not clusters_differ :
        anc_region = range(int(nSNPs *(1-prop_causal)))
        causal_region = range(max(anc_region)+1, nSNPs)

        anc_regionsize = sum(maf_filter[anc_region])
        causal_regionsize = sum(maf_filter[causal_region])
        causal_snps = [False for i in range(anc_regionsize)]
        causal_snps += [True for i in range(causal_regionsize)]

    # Use MAF filter
    filtered_geno = genotypes[:,  maf_filter]
    
    # standardize the genotpyes
    allele_freqs = np.mean(filtered_geno, axis=0) / 2
    filtered_geno = np.matrix((filtered_geno- 2 * allele_freqs) /
                          np.sqrt(2 * allele_freqs * (1 - allele_freqs)))

    # Calc standardized GRM
    GRM = np.dot(filtered_geno, filtered_geno.T) / nSNPs
    
    
    
    # project GRM onto pc space
    pcs = pd.DataFrame(PCA(n_components = 20).fit_transform(np.asarray(GRM)))
    pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
    
    return genotypes, GRM, pcs, causal_snps


