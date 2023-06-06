#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 13:08:24 2023

@author: christian
"""
import numpy as np
import pandas as pd

def sim_pop_alleles(seed, theta_alleles = [0.8, 0.2], nclusts=1, nSNPs = 1000, shared_causal= 0.8, shared_noncausal = 0.8, prop_causal = 0.1, maf_filter = 0.1) :
    """
    Simulate the allele frequencies for a common ancestor and for all genetic clusters taking into account if the causal snps are shared or not.

    Parameters
    ----------
    theta_alleles : list of floats, optional
        between 0, 1, parameter controling the similarity between genetic clusters. near 0- similar clusters, near 1- distinct clusters. 
        First value is for nonshared regions, second is for shared regions
        The default is 0.5.
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
    if theta_alleles[0] < theta_alleles[1] :
        raise ValueError("Shared regions should be more conserved than nonshared regions (i.e. theta_alleles[1] should be less than theta_alleles[0]")
    
    rng = np.random.default_rng(seed)
    boundaries = {"causal_shared" : np.array([0, int(prop_causal * nSNPs * shared_causal)])}
    boundaries["causal_nonshared"] = np.array([0 , int(prop_causal * nSNPs * (1-shared_causal))]) + np.max(boundaries["causal_shared"]) 
    boundaries["noncausal_shared"] = np.array([0, int((1-prop_causal) * nSNPs * shared_noncausal)]) + np.max(boundaries["causal_nonshared"])
    boundaries["noncausal_nonshared"] = np.array([0, int((1-prop_causal) * nSNPs * (1-shared_noncausal))]) + np.max(boundaries["noncausal_shared"])
    
    
    # simulate ancestral frequency of each SNP, dimensions = (SNPs,)
    ancest_freqs = rng.uniform(low=0.1, high=0.9, size=nSNPs)

    if nclusts == 1:
        cluster_freqs = ancest_freqs

    else :

        # simulate allele frequencies for each cluster
        cluster_freqs = np.zeros((nclusts, nSNPs))
        
        # simulate the allele frequencies across the four regions and join them for each cluster
        for (region, bounds) in boundaries.items() :
            bounds = range(*bounds)
            if "nonshared" in region :
                theta = theta_alleles[0]
            else :
                theta = theta_alleles[1]
    
            cluster_freqs[:, bounds] = rng.beta(ancest_freqs[bounds] * (1- theta) / theta,
                                      (1-ancest_freqs[bounds]) *
                                      (1-theta)/theta,
                                      size=(nclusts, len(bounds)))
    
    # Get vectors annotating genomes for being shared or causal
    shared_idx = list(range(*boundaries["causal_shared"])) + list(range(*boundaries["noncausal_shared"]))
    causal_idx = list(range(*boundaries["causal_shared"])) + list(range(*boundaries["causal_nonshared"]))
    
    if nclusts >1 :
        # filter out the causal_idx for which the cluster_freqs are less than the maf or greater than 1-maf for at least one columns
        causal_idx = np.array(causal_idx)[np.all((cluster_freqs[:, causal_idx] > maf_filter) & (cluster_freqs[:, causal_idx] < 1-maf_filter), axis=0)]
    else :
        causal_idx = np.array(causal_idx)[(cluster_freqs[causal_idx] > maf_filter) & (cluster_freqs[causal_idx] < 1-maf_filter)]

    return ancest_freqs, cluster_freqs, np.array(shared_idx), np.array(causal_idx)

def assign_clusters(df, rng, theta_alleles=0.5, nclusts=1, nsites = 1, site_comp="IID", dominance=2, eq_sites = False):
    """
    Simulate the genetic cluster that each subject will be observed from 

    Parameters
    ----------
    df : pandas datafrmae containing site assignments
        containing site assignments in one column.
    theta_alleles : float, optional
        between 0, 1, parameter controling the similarity between genetic clusters. near 0- similar clusters, near 1- distinct clusters.
        The default is 0.5.
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
    # site_comp = ["EQUAL", "RAND", "IID", "HET"]
    if eq_sites == True:
        site_comp = "EQUAL"
    # Theta alleles controls how much ancstry is conserved. low theta_alleles = high conservation
    theta_alleles = theta_alleles
    # Clusts are distinct ancestries
    nclusts = nclusts
    site_comp = site_comp


    # Sample ancesrties for each individual Dim = (nsubjects,)
    if site_comp == "EQUAL":
        # Make sure each cluster is equally represented i neach site
        subj_ancestries = np.tile(np.arange(nclusts), int(
            nsubjects  / nclusts))

    elif site_comp == "RAND":
        subj_ancestries = rng.integers(
            low=0, high=nclusts, size=nsubjects)

    elif site_comp == "IID":
        # Sample probabilities of ancestries for each site (nclusters x nsites) from same dirichlet distribution
        pop_probs = rng.dirichlet(
            np.repeat(1, nclusts), size= nsites)
        subj_ancestries = []
        for s in df["abcd_site"]:
            # Sample subject i's ancestry given the proportion at site s
            subj_ancestries.append(
                rng.choice(np.arange(nclusts), p=pop_probs[s]))

    elif site_comp == "HET":
        pop_probs = np.zeros((nsites, nclusts))
        # sample each site largely from  one cluster by giving large alpha to one cluster for each site
        alphas = np.ones((nsites, nclusts))
        for i in range(alphas.shape[0]):
            # which cluster will be high proportion
            high_prop = rng.choice(nclusts)
            alphas[i, high_prop] = dominance
            pop_probs[i] = rng.dirichlet(alphas[i])

            subj_ancestries = []
        for s in df["abcd_site"]:
            # Sample subject i's ancestry given the proportion at site s
            subj_ancestries.append(
                rng.choice(np.arange(nclusts), p=pop_probs[s]))

    return subj_ancestries

    # Checked samples came from ancester with the following
    # plt.scatter(sim1.ancest_freqs, sim1.pop_freqs.mean(axis = 0))
