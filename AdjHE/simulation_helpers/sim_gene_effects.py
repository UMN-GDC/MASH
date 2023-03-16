#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:42:57 2023

@author: christian
"""
import numpy as np


def sim_gen_effects(rng, genotypes, causals, df, prop_causal=0.1, site_dep=False, clusters_differ= False):
    """
    

    Parameters
    ----------
    rng : TYPE
        DESCRIPTION.
    genotypes : TYPE
        DESCRIPTION.
    causals : TYPE
        DESCRIPTION.
    df : TYPE
        DESCRIPTION.
    prop_causal : TYPE, optional
        DESCRIPTION. The default is 0.1.
    site_dep : TYPE, optional
        DESCRIPTION. The default is False.
    clusters_differ : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    Gene_contrib : TYPE
        DESCRIPTION.

    """
    nsubjects, nSNPs = genotypes.shape
    nCausal = int(nSNPs * prop_causal)
    nsites = len(np.unique(df["abcd_site"]))
    # select causal snos
    causals = rng.choice(nSNPs, nCausal, replace=False, shuffle=False)

    # Select just causal genes
    if not clusters_differ :
        Xcausal = np.matrix(genotypes[:, causals])
        nCausal = len(causals)
    elif clusters_differ :
        Xcausal = np.matrix(genotypes[:, causals])

    if site_dep:
        Gene_contrib = np.zeros((nsubjects, ))
        # Simulate Causal effects for SNP's based on site
        causal_eff = rng.normal(0, 1, size=(nCausal, nsites))

        for i in range(nsubjects):
            # determine the persons site
            site = df["abcd_site"][i]

            # multiply each persons genotype with the snp effects for the specified site
            Gene_contrib[i, ] = (
                np.dot(Xcausal[i, :], causal_eff[:, site]))

    else:
        # sim effect from each SNP (Sample from N(0,1), later rescale to get desired variance contributions)
        causal_eff = rng.normal(0, 1, (Xcausal.shape[1], 1))
        Gene_contrib = np.array(Xcausal * causal_eff).flatten()

    return Gene_contrib
