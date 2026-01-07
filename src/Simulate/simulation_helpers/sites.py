#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
import numpy as np
import pandas as pd


def sim_sites(rng, nsubjects = 1000, nsites=1, siteDistribution="EQUAL", random_BS = True, nphenos = 1):
    """
    Simulate a site vector for the specified number of subjects and sites. Also, control wether the sites will be of equal sites or 
    sampled.

    Parameters
    ----------
    rng : numpy.random._generator.Generator 
        Random number generator to increase reproducibility
    nsubjects : int, optional
        Number of subjects to simulate. The default is 1000.
    nsites : int, optional
        Number of sites to simulate. The default is 1.
    eq_sites : bool, optional
        True- make all sites equal sized (note that the number of subjects must be divisible by the number of sites). The default is False.
    random_BS : bool, optional
        Whether the effecs (betas) for each site is randomly generated or not. The default is True.

    Returns
    -------
    a pandas dataframe (nsubjects x 2) numpy array with integers specifyin the simulated site ID and the associated effect of the site on
    the phenotype

    """
    # Site assignment
    if siteDistribution == "EQUAL":
        # For equal assignment (n needs to be divisible by nsites * nclusts)
        Groups = np.repeat(
            np.arange(start=1, stop=nsites + 1), int(nsubjects/nsites))
    else:
        # For unequal assignment sample randomly from uniform
        Groups = rng.choice(nsites, size= nsubjects, replace=True)

        # Simulate site effects (will rescale to desired contribution later)
    if random_BS :
        Bs = np.matrix(np.random.normal(0,  1, (nsites, nphenos)))
    else :
        # make site effects evenly spaced in a nsites x nphenos matrix 
        Bs = np.arange(nsites * nphenos, dtype = "float64", like = np.zeros(nsites, nphenos) )
        Bs -= Bs.mean()
        Bs=  np.reshape(Bs, (nsites, nphenos))
    if nsites == 1:
        Bs = np.array([0])
        site_contribs = 0
    else : 
        # Make a dummy matrix
        Groups_dum = np.matrix(pd.get_dummies(Groups))
    
        # dot product between Groups_dum and Bs save them as a dataframe with columns 
        # Site_contrib followed by the number of the site
        site_contribs = np.dot(Groups_dum, Bs)
        site_contribs = pd.DataFrame(site_contribs, columns = [f"Site_contrib{i}" for i in range(nphenos)])
    
        # Add the groups as a column to site_contribs dataframe
        site_contribs["abcd_site"] = Groups
        # make abcd_site the first column
        site_contribs = site_contribs[["abcd_site"] + [col for col in site_contribs.columns if col != "abcd_site"]]
    return site_contribs 

