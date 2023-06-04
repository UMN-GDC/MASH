#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
import logging
import numpy as np
import pandas as pd
from pandas_plink import read_plink
from Simulator.simulation_helpers.sites import sim_sites
from Simulator.simulation_helpers.clusters import sim_pop_alleles, assign_clusters
from Simulator.simulation_helpers.genos import sim_genos
from Simulator.simulation_helpers.gene_effects import sim_gen_effects
from Simulator.simulation_helpers.covariates import sim_covariates
from Simulator.simulation_helpers.pheno import sim_pheno
from Simulator.simulation_helpers.plink_pheno import sim_plink_pheno


class pheno_simulator():
    def __init__(self, nsubjects=1000, nSNPs=1000, plink_prefix=None):
        self.plink_prefix = plink_prefix
        
        if plink_prefix == None :
            self.nsubjects = nsubjects
            self.nSNPs = nSNPs

            # Seed the simulated dataframe
            self.df = pd.DataFrame({"FID": np.arange(start=1, stop=nsubjects + 1), 
                                    "IID": np.arange(start=1, stop=nsubjects + 1)})

            self.rng = np.random.default_rng()
        else :
            self.plink_prefix= plink_prefix
            
            (bim, fam, bed) = read_plink(plink_prefix)
            (self.nSNPs, self.nsubjects) = bed.shape
            self.df = pd.DataFrame({"FID": fam.fid,
                                    "IID": fam.iid})
            # bed is (nSNPs x nsubjects) so transpose it
            self.genotypes = bed.T
            self.rng = np.random.default_rng()


    def sim_sites(self, nsites=1, eq_sites=False, random_BS = True):
        self.eq_sites = eq_sites
        self.nsites = nsites
        self.df[["abcd_site", "Site_contrib"]] = sim_sites(rng = self.rng, nsubjects = self.nsubjects, 
                                                           nsites=nsites, eq_sites=eq_sites, random_BS = random_BS)

    def sim_pops(self, theta_alleles=[0.8, 0.2], nclusts=1, site_comp="IID", dominance=2, 
                 shared_causal= 0.8, shared_noncausal = 0.8, prop_causal = 0.1):
        # site_comp = ["EQUAL", "RAND", "IID", "HET"]
        if self.eq_sites == True:
            site_comp = "EQUAL"
        # Theta alleles controls how much ancstry is conserved. low theta_alleles = high conservation
        self.theta_alleles = theta_alleles
        # Clusts are distinct ancestries
        self.nclusts = nclusts
        self.site_comp = site_comp
        self.shared_causal= shared_causal
        self.shared_noncausal = shared_noncausal
        
        # dominance is only referenced for heterogeneous cluster sampling
        if site_comp == "HET" :
            self.dominance = dominance

        # Simualte allele frequencies for common ancestor and for genetic clusters
        (self.ancest_freqs, self.cluster_frequencies, self.shared_idx, self.causal_idx) = sim_pop_alleles(seed = self.rng, theta_alleles = self.theta_alleles,
                                                                                                          nclusts=self.nclusts, nSNPs = self.nSNPs,
                                                                                                          shared_causal= shared_causal, shared_noncausal = shared_noncausal,
                                                                                                          prop_causal = prop_causal)

        # Sample ancesrties for each individual Dim = (nsubjects,)
        self.df["subj_ancestries"] = assign_clusters(df = self.df, rng = self.rng, theta_alleles=self.theta_alleles, nclusts=self.nclusts,
                                                     nsites = self.nsites, site_comp= site_comp, dominance=dominance, eq_sites = self.eq_sites)

    def sim_genos(self, prop_causal=0.1, maf_filter = 0.05, admixing = False):
        self.prop_causal = prop_causal
        self.maf_filter  = maf_filter
        self.admixing = admixing
        (self.genotypes, self.GRM, pcs) = sim_genos(seed = self.rng, ancestral_frequencies= self.ancest_freqs,
                                                    cluster_frequencies = self.cluster_frequencies, 
                                                    subject_ancestries = self.df["subj_ancestries"],
                                                    prop_causal=self.prop_causal, admixing = admixing)
        
        # update the number of SNPs
        self.nSNPS_post_filter = self.genotypes.shape[1]
        self.df = pd.concat([self.df, pcs], axis=1)


    def sim_gen_effects(self, clusters_differ = False, variance_propto_frequency = True, maf_filter =0.1):
        # genetic contribution
        self.clusters_differ = clusters_differ
        self.variance_propto_frequency = variance_propto_frequency
        self.maf_filter = maf_filter
        # simulate sared genetic contribution
        self.df["Gene_contrib"], self.causals = sim_gen_effects(rng = self.rng, genotypes = self.genotypes, causals= self.causal_idx,
                                                  prop_causal=self.prop_causal, variance_propto_frequency = variance_propto_frequency,
                                                  maf_filter = maf_filter) 

        # if differ is true then simulate clusters specific genetic contribution
        #if clusters_differ :
        if False :
            for cluster in range(self.nclusts) :
                # find the indices of the group variable
                idx = self.df["subj_ancestries"] == cluster

                # subset the genotype rows by the above idx
                temp_genos = self.genotypes[idx, :]
                # Simulate the contribution from eahc cluster 
                temp = np.zeros(self.df.shape[0])

                temp[self.df.subj_ancestries == cluster] = sim_gen_effects(rng=self.rng, genotypes = temp_genos, causals= self.causal_idx,
                                                                             prop_causal=self.prop_causal, variance_propto_frequency = variance_propto_frequency)
                self.df["Gene_contrib_c_" + str(cluster)] = temp



    def sim_covars(self, cov_effect= True, ortho_cov = False) :
        self.cov_effect= cov_effect
        self.ortho_cov = ortho_cov
        
        self.df["Xc"], self.df["Covar_contrib"] = sim_covariates(rng = self.rng, nclusts=self.nclusts, 
                                                                 df = self.df, cov_effect=cov_effect,ortho_cov = ortho_cov)
        
    def sim_pheno(self, var_comps=[0.5, 0.25, 0.25], phen = 1, site_het = False, cluster_contribs = None):
        self.cluster_contribs=  cluster_contribs

        self.var_comps = var_comps
        self.site_het = site_het
        pheno_contribs  = sim_pheno(rng = self.rng, df = self.df, var_comps=var_comps, 
                                    phen = phen, site_het = site_het, nsites = self.nsites,nclusts =self.nclusts,
                                    cluster_contribs = cluster_contribs)
        
        # join pheno_contribs dataframe to simulated dataframe as new columns
        self.df.update(pheno_contribs)
            
    def full_sim(self, sigma = [0.5,0.25,0.25], site_comp="IID",
                 nsites=30, theta_alleles=[0.5, 0.5], nclusts=5, dominance=5,
                 prop_causal=0.25, nsubjects=1000,
                 nnpc=0, phens = 2, site_het = False, cov_effect = True,
                 ortho_cov = True, random_BS = True, npcs = 0, admixing = False, maf_filter= 0.1, clusters_differ = False):
        
        # Only simulate genes if plink file not specified
        if self.plink_prefix == None: 
  
            if site_comp == "EQUAL" :
                eq_sites = True
            else :
                eq_sites= False
            
            # Run through full simulation and estimation
            self.sim_sites(nsites= nsites, eq_sites=eq_sites, random_BS = random_BS)
            
            self.sim_pops(theta_alleles=theta_alleles, nclusts=nclusts, site_comp= site_comp, dominance=dominance)
            self.sim_genos(prop_causal=prop_causal, admixing = admixing)
            self.sim_gen_effects(maf_filter = maf_filter, clusters_differ = clusters_differ)
            self.sim_covars(cov_effect= cov_effect, ortho_cov = ortho_cov)
            
            for i in range(phens) :
                self.sim_pheno(var_comps=sigma, phen = i, site_het = site_het)
        
        else : 
            for i in range(phens) :
                phenoname = "Y" + str(i)
                self.df[phenoname] = sim_plink_pheno(rng = self.rng, bed = self.genotypes, sigma= sigma, prop_causal = prop_causal, npcs=npcs)

    def summary(self) :
        """
        Summarize all of the input parameters attributes of the simulation object

        Returns
        -------
        None.

        """
        logging.info("Simulated data summary")
        logging.info("Number of subjects: ", self.nsubjects)
        logging.info("Number of SNPs: ", self.nSNPs)
        logging.info("Number of clusters: ", self.nclusts)
        logging.info("Proportion of causal SNPs: ", self.prop_causal)
