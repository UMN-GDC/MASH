#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 13:08:24 2023

@author: christian
"""
import numpy as np
import pandas as pd

def sim_pop_alleles(rng, theta_alleles = [0.8, 0.2], nclusts=1, nSNPs = 1000, shared_causal= 0.8, shared_noncausal = 0.8, prop_causal = 0.1, maf_filter = 0.1) :
    """
    Simulate the allele frequencies for a common ancestor and for all genetic clusters taking into account if the causal snps are shared or not.

    Parameters
    ----------
    theta_alleles : list of floats, optional
        between 0, 1, parameter controling the similarity between genetic clusters. near 0- similar clusters, near 1- distinct clusters. 
        First value is for nonshared regions, second is for shared regions
        The default is [0.8, 0.2].
    nclusts : int, optional
        number of genetic clusters. The default is 1.
    nSNPs : int, optional
        Number of snps to simulate for the genome. The default is 1000.
    shared_causal : floa64
        specify the proportion of causal snps shared between genetic clusters, between 0 and 1.
    shared_noncasaul : float64
        specify the proportion of non-causal snps shared between genetic clusters, between 0 and 1.
    prop_causal : float64
        specify the proportion snps that are causal.
    maf_filter : float64
        between 0 and 0.5, the minimum allowable maf for causal snps

    Returns
    -------
    tuple with (nSNPs x 1) numpy array where the first column corresponds to the common ancestral frequencies,
    (nSNPs x (nclusts +1)) numpy array columns correspond to each subclusters allele frequencies.
    indices of shared snps
    incides of causal snps

    """
    rng = np.random.default_rng(rng)
    theta = theta_alleles
    # simulate ancestral frequency of each SNP, dimensions = (SNPs,)
    ancest_freqs = rng.uniform(low=0.1, high=0.9, size=nSNPs)
    cluster_freqs = rng.beta(ancest_freqs * (1- theta) / theta,
                                      (1-ancest_freqs) *
                                      (1-theta)/theta,
                                      size=(nclusts, nSNPs))
    return ancest_freqs, cluster_freqs 


def assign_clusters(df, rng, nclusts=1, nsites = 1):
    """
    Simulate the genetic cluster that each subject will be observed from 

    Parameters
    ----------
    df : pandas datafrmae containing site assignments
        containing site assignments in one column.
    nclusts : int, optional
        number of genetic clusters. The default is 1.
    nsites : int, optional
        number of sites. The default is 1.
    site_comp : str, optional
        Specifying the balance of genetic clusters within each site. One of ["IID", "EQUAL", "RAND", "HET"]. The default is "IID".
    dominance : float, optional
        if site_comp is HET, this specifies the change in the sampling exponent associated with the most prevelant genetic cluster.
        The default is 2.
    eq_sites : bool, optional
        If true, forces all sites to have equal balance with all genetic clusters. The default is False.

    Returns
    -------
    (nsubjects x 1) numpy array specifying the cluster that each subject will be sampled from

    """
    nsubjects= df.shape[0]
    # Make sure each cluster is equally represented i neach site
    subj_ancestries = np.tile(np.arange(nclusts), int(
        nsubjects  / nclusts))

    return subj_ancestries

    # Checked samples came from ancester with the following
    # plt.scatter(sim1.ancest_freqs, sim1.pop_freqs.mean(axis = 0))
