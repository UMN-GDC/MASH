#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:09:08 2023

@author: christian
"""
import numpy as np
import statsmodels.formula.api as smf
from Simulator.simulation_helpers.sim_effects import sim_effects 

def sim_pheno(rng, df, var_comps=[0.5, 0.25, 0.25], phen = 1, site_het = False, nsites = 1, nclusts =1, cluster_contribs = None):
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
    # make sure no zero variance components
    for i, v in enumerate(var_comps):
        if v == 0:
            var_comps[i] = 1e-6

    # Sim errors
    if site_het :
        # Sample the site variances
        site_var = rng.gamma(4, 4, nsites)
        # Sample error from the specified variance
        errors = []
        for i in range(nsubjects):
            # determine the persons site
            site = df["abcd_site"][i]
            # Sample their errors given that sites' variance
            errors += [rng.normal(0, site_var[site -1 ])]
    else : 
        errors = rng.normal(0, 1, size=nsubjects)

    # Calculate desired variance ratios
    StoG = var_comps[1]/var_comps[0]
    EtoG = var_comps[2]/var_comps[0]
    
    if nclusts > 1 : 
        # Get pc column names
        pcs = ["pc_" + str(i + 1) for i in range(nclusts)]
        # Build regression equation
        form= "Gene_contrib ~ " + " + ".join(pcs)
        # Find the Genetic_contribution after accountring for race
        resid = smf.ols(formula = form, data = df).fit().resid
        # Find the emprical variance after accounting for race
        gen_var = np.var(resid)
    else :
        gen_var = np.var(df["Gene_contrib"])
        
    # Calculate empricial variance due to site
    site_var = np.var(df["Site_contrib"])

    # Find the empirical ratio of variances
    StoG_sim = site_var / gen_var
    site_variance_scaling = StoG_sim / StoG

    EtoG_sim = np.var(errors) / gen_var
    error_variance_scaling = EtoG_sim / EtoG

    if cluster_contribs != None :
        for (cluster, cluster_contrib) in enumerate(cluster_contribs) :
            CtoG = cluster_contrib / var_comps[0]
            CtoG_sim = np.var(df["Gene_contrib_c" + str(cluster)])
            Cvariance_scaling = CtoG_sim/CtoG
            df["Gene_contrib_c" + str(cluster)] = df["Gene_contrib_c" + str(cluster)] / np.sqrt(Cvariance_scaling)

        return_columns = ["Gene_contrib", "Site_contrib", "errors"] + ["Gene_contrib_c" + str(i) for i in range(len(cluster_contribs))]
    else : 
        return_columns = ["Gene_contrib", "Site_contrib", "errors"] 


    # Scale site effects so total contributed variance is as presecribed
    df["Site_contrib"] = rescalar(df["Gene_contrib"], df["Site_contrib"], var_comps[0], var_comps[1])
    df["Site_contrib"] = df["Site_contrib"] / np.sqrt(site_variance_scaling)
    site_var = np.var(df["Site_contrib"])

    # Sample errors from normal scaled by the ratio of the intedended variance between genetic and error effects
    df["errors"] = (errors / np.sqrt(error_variance_scaling)).flatten()
    # Sim phenotype
    if phen==0 :
        phenoname = "Y"
    else :
        phenoname = "Y" + str(phen)

    if np.sum(np.isnan(df.Site_contrib)) > 100 :
        df.Site_contrib = 0
    df[phenoname] = df[["Gene_contrib", "Site_contrib", "errors", "Covar_contrib"]].sum(axis=1)
    df[phenoname] = df[phenoname] - np.nanmean(df[phenoname])
    
    return df[return_columns + [phenoname]]
