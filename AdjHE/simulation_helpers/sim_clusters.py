#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 13:08:24 2023

@author: christian
"""
import numpy as np
import pandas as pd

def sim_pop_alleles(rng, theta_alleles = 0.5, nclusts=1, nSNPs = 1000) :
    """
    Simulate the allele frequencies for a common ancestor and for all genetic clusters

    Parameters
    ----------
    theta_alleles : float, optional
        between 0, 1, parameter controling the similarity between genetic clusters. near 0- similar clusters, near 1- distinct clusters.
        The default is 0.5.
    nclusts : int, optional
        number of genetic clusters. The default is 1.
    nSNPs : int, optional
        Number of snps to simulate for the genome. The default is 1000.

    Returns
    -------
    tuple with (nSNPs x 1) numpy array where the first column corresponds to the common ancestral frequencies,
    (nSNPs x (nclusts +1)) numpy array columns correspond to each subclusters allele frequencies.

    """
    # simulate ancestral frequency of each SNP, dimensions = (SNPs,)
    ancest_freqs = rng.uniform(low=0.1, high=0.9, size=nSNPs)

    # simulate allele frequency for each population, Dim = (nclusts x SNPS)
    cluster_frequencies = rng.beta(ancest_freqs * (1- theta_alleles) / theta_alleles,
                              (1-ancest_freqs) *
                              (1-theta_alleles)/theta_alleles,
                              size=(nclusts, nSNPs))
    
    return ancest_freqs, cluster_frequencies



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
