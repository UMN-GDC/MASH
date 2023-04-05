#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 13:29:28 2023

@author: christian
"""

import numpy as np

def sample_admixed_genotype(seed, cluster_freqs, admixture, nsubjects=1) :
    """
    Sample subject genotypes from admixtures of multiple genetic clusters

    Parameters
    ----------
    cluster_freqs : numpy array
        (#Clusters x #SNPs) array speciying the frequency for each SNP for each genetic cluster.
    nsubjects : int
        the number of subjects to simulate with the given admixing.
    admixture : numpy array
        1-D array specying the propoertion of allele frequencies coming from each population.

    Returns
    -------
    (#subjects x #SNPs) array with subject genotypes.

    """
    rng = np.random.default_rng(seed)
    (nclusts, nSNPs) = cluster_freqs.shape
    # Seed empty probability vector
    admix_freq = np.zeros(nSNPs)
    
    # Sample one population frequency from each cluster
    for snp in range(nSNPs) :
        # Define allele frequency for unique admixture of given subject
        admix_freq[snp] = rng.choice(cluster_freqs[:,snp], p = admixture)
        
    # return that subjects phenotype
    return rng.binomial(np.repeat(2, nSNPs), admix_freq, size = (nsubjects, nSNPs))

