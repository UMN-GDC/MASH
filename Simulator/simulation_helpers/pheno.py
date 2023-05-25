#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:09:08 2023

@author: christian
"""
import numpy as np
import statsmodels.formula.api as smf
from Simulator.simulation_helpers.sim_effects import sim_effects 

# Function: Takes vector X, and scales it such that its variance is equal to V
def rescalar(X, V): 
    """
    Rescale the vector X such that its variance is equal to V.

    Parameters
    ----------
    X : numpy array
        vector to be rescaled.
    V : float
        desired variance.

    Returns
    -------
    numpy array
        rescaled vector.

    """
    return X * np.sqrt(V / np.nanvar(X))


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


    if nclusts > 1 :
        # If clusters exist, separate PC effects from the local genetic effects
        # Get pc column names
        pcs = ["pc_" + str(i + 1) for i in range(nclusts)]
        # Build regression equation
        form= "Gene_contrib ~ " + " + ".join(pcs)
        # Find the Genetic_contribution after accountring for race
        mod = smf.ols(formula = form , data= df).fit()
        df["Gene_contrib"] = mod.resid
        df["PC_contrib"] = mod.predict()
    else :
        # If no clusters, there is no PC contribution 
        df["PC_contrib"] = 0

        
    if cluster_contribs != None :
        for (cluster, cluster_contrib) in enumerate(cluster_contribs) :
            CtoG = cluster_contrib / var_comps[0]
            CtoG_sim = np.var(df["Gene_contrib_c" + str(cluster)])
            Cvariance_scaling = CtoG_sim/CtoG
            df["Gene_contrib_c" + str(cluster)] = df["Gene_contrib_c" + str(cluster)] / np.sqrt(Cvariance_scaling)

        return_columns = ["Gene_contrib", "PC_contrib", "Site_contrib", "errors"] + ["Gene_contrib_c" + str(i) for i in range(len(cluster_contribs))]
    else : 
        return_columns = ["Gene_contrib", "PC_contrib", "Site_contrib", "errors"] 


    # Scale genetic, site and error contributions to get the desired heritability
    df["Gene_contrib"] = rescalar(df["Gene_contrib"], var_comps[0])
    df["Site_contrib"] = rescalar(df["Site_contrib"], var_comps[1])
    df["errors"] = rescalar(errors, var_comps[2]) 
    
    if phen == 0:
        phenoname = "Y"
    else : 
        phenoname = "Y" + str(phen-1)

    if np.sum(np.isnan(df.Site_contrib)) > 10 :
        df.Site_contrib = 0
    df[phenoname] = df[["Gene_contrib", "PC_contrib", "Site_contrib", "errors", "Covar_contrib"]].sum(axis=1)
    df[phenoname] = df[phenoname] - np.nanmean(df[phenoname])
    
    return df[return_columns + [phenoname]]
