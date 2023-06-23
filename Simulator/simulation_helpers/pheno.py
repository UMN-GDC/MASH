#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:09:08 2023

@author: christian
"""
import numpy as np
import statsmodels.formula.api as smf
from Simulator.simulation_helpers.sim_effects import sim_effects 



def sim_pheno(rng, df, h2= 0.5, phen = 1, nsites = 1, nclusts =1):
    """
    Simulate the phenotype given the differnt contributions from site, genetic, and error, and scale them by the var_comps variable.

    Parameters
    ----------
    rng : numpy random number generator
        numpy random number generator.
    df : pandas dataframe
        dataframe containing contributions from genes and sites.
    var_comps : list of floats, optional
        specify the proportion of variance attributable to genetics, site, and error, respectively. The default is [0.5, 0.25, 0.25].
    phen : int, optional
        specify how many phenotypes to simulate. The default is 1.
    site_het : bool, optional
        should the sites have heterogeneous error. The default is False.
    nsites : int, optional
        number of sites in the simulation. The default is 1.
    nclusts : int, optional
        number of clusters in thhe study. The default is 1.
    cluster_contribs : list of floats, optional
        A list specifying additional cluster specific effects. If these values differ, they will change the heritability of each 
        cluster (they will be added to the shared variance component in var_comps)

    Returns
    -------
    pandas series
        simulated genetic contribution as a pandas series.
    pandas series
        simulated site effect as a pandas series.
    pandas series
        simulated errors as a pandas series.
    pandas series
        simulated phenotype as a pandas series.
    """
    nsubjects = df.shape[0]
    df["VG"] = df.groupby("subj_ancestries")["Gene_contrib"].transform(np.nanvar)
    df["errors"] = rng.normal(np.repeat(0, nsubjects), np.sqrt(df["VG"] * (1- h2)/h2)) 
   
    if phen == 0:
        phenoname = "Y"
    else : 
        phenoname = "Y" + str(phen-1)

    if np.sum(np.isnan(df.Site_contrib)) > 2 :
        df.Site_contrib = 0
    df[phenoname] = df[["Gene_contrib", "Site_contrib", "errors", "Covar_contrib"]].sum(axis=1)
    df[phenoname] = df[phenoname] - np.nanmean(df[phenoname])
    return df[["Gene_contrib", "Site_contrib", "errors", "Covar_contrib", phenoname]]
