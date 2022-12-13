#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
#from plotly.offline import plot
#import plotly.express as px
import os
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


rng = np.random.default_rng(123)


# TO DO:
# Compare with naively splitting into sites and estimating (see if heritability estimate changes by site)


# %%
class pheno_simulator():
    def __init__(self, nsubjects, nSNPs):
        self.nsubjects = nsubjects
        self.nSNPs = nSNPs

        # Seed the dataframe
        self.df = pd.DataFrame({"FID": np.arange(start=1, stop=nsubjects + 1),
                                "IID": np.arange(start=1, stop=nsubjects + 1),
                                })

    def sim_sites(self, nsites=1, eq_sites=False):
        self.eq_sites = eq_sites
        # Site assignment
        if eq_sites:
            # For equal assignment (n needs to be divisible by nsites * nclusts)
            Groups = np.repeat(
                np.arange(start=1, stop=nsites + 1), int(self.nsubjects/nsites))
        else:
            # For unequal assignment sample randomly from uniform
            Groups = rng.choice(nsites, size=self.nsubjects, replace=True)

        self.df["abcd_site"] = Groups
        # Simulate site effects (will rescale to desired contribution later)
        Bs = np.matrix(np.random.normal(0,  1, nsites)).T
        # Make a dummy matrix
        Groups_dum = np.matrix(pd.get_dummies(self.df["abcd_site"]))
        self.df["Site_contrib"] = np.array(Groups_dum * Bs).flatten()
        self.nsites = nsites

    def sim_pops(self, theta_alleles=0.5, nclusts=1, site_comp="IID", dominance=2):
        # site_comp = ["EQUAL", "RAND", "IID", "HET"]
        if self.eq_sites == True:
            site_comp = "EQUAL"
        # Theta alleles controls how much ancstry is conserved. low theta_alleles = high conservation
        self.theta_alleles = theta_alleles
        # Clusts are distinct ancestries
        self.nclusts = nclusts
        self.site_comp = site_comp

        # simulate ancestral frequency of each SNP, dimensions = (SNPs,)
        self.ancest_freqs = rng.uniform(low=0.1, high=0.9, size=self.nSNPs)

        # simulate allele frequency for each population, Dim = (nclusts x SNPS)
        self.pop_freqs = rng.beta(self.ancest_freqs * (1-self.theta_alleles) / self.theta_alleles,
                                  (1-self.ancest_freqs) *
                                  (1-self.theta_alleles)/self.theta_alleles,
                                  size=(self.nclusts, self.nSNPs))

        # Sample ancesrties for each individual Dim = (nsubjects,)
        if site_comp == "EQUAL":
            # Make sure each cluster is equally represented i neach site
            subj_ancestries = np.tile(np.arange(self.nclusts), int(
                self.nsubjects  / self.nclusts))

        elif site_comp == "RAND":
            subj_ancestries = rng.integers(
                low=0, high=self.nclusts, size=self.nsubjects)

        elif site_comp == "IID":
            # Sample probabilities of ancestries for each site (nclusters x nsites) from same dirichlet distribution
            pop_probs = rng.dirichlet(
                np.repeat(1, self.nclusts), size=self.nsites)
            subj_ancestries = []
            for s in self.df["abcd_site"]:
                # Sample subject i's ancestry given the proportion at site s
                subj_ancestries.append(
                    rng.choice(np.arange(self.nclusts), p=pop_probs[s]))

        elif site_comp == "HET":
            pop_probs = np.zeros((self.nsites, self.nclusts))
            # sample each site largely from  one cluster by giving large alpha to one cluster for each site
            alphas = np.ones((self.nsites, self.nclusts))
            for i in range(alphas.shape[0]):
                # which cluster will be high proportion
                high_prop = rng.choice(self.nclusts)
                alphas[i, high_prop] = dominance
                pop_probs[i] = rng.dirichlet(alphas[i])

                subj_ancestries = []
            for s in self.df["abcd_site"]:
                # Sample subject i's ancestry given the proportion at site s
                subj_ancestries.append(
                    rng.choice(np.arange(self.nclusts), p=pop_probs[s]))

        self.df["subj_ancestries"] = subj_ancestries

        # Checked samples came from ancester with the following
        # plt.scatter(sim1.ancest_freqs, sim1.pop_freqs.mean(axis = 0))

    def sim_genos(self):
        if self.nclusts == 1:
            # simulate genotypes
            genotypes = rng.binomial(n=np.repeat(2, self.nSNPs), p=self.pop_freqs[0],
                                     size=(self.nsubjects, self.nSNPs))
        else:
            # simulate genotypes
            genotypes = rng.binomial(n=np.repeat(2, self.nSNPs), p=self.pop_freqs[self.df["subj_ancestries"]],
                                     size=(self.nsubjects, self.nSNPs))
        # keep SNPs with MAF greater than 0.05
        maf_filter = np.logical_and((np.sum(genotypes, axis=0) / (2 * self.nsubjects)) > 0.05,
                                    (np.sum(genotypes, axis=0) / (2 * self.nsubjects)) < 0.95)
        self.genotypes = genotypes[:,  maf_filter]

        # get new number of SNPs
        self.nSNPs = self.genotypes.shape[1]

        # reference allele frequ
        # standardize the genotpyes
        allele_freqs = np.mean(self.genotypes, axis=0) / 2
        genotypes = np.matrix((self.genotypes - 2 * allele_freqs) /
                              np.sqrt(2 * allele_freqs * (1 - allele_freqs)))
        self.GRM = np.dot(genotypes, genotypes.T) / self.nSNPs
        # project GRM onto pc space
        pcs = pd.DataFrame(PCA(n_components=20).fit_transform(self.GRM))
        pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
        # add the pcs to the dataframe
        self.df = pd.concat([self.df, pcs], axis=1)


    def sim_gen_effects(self, prop_causal=0.1, site_dep=False):
        nCausal = int(self.nSNPs * prop_causal)
        # select causal snos
        causals = rng.choice(self.nSNPs, nCausal, replace=False, shuffle=False)

        # Select just causal genes
        Xcausal = np.matrix(self.genotypes[:, causals])

        if site_dep:
            Gene_contrib = np.zeros((self.nsubjects, ))
            # Simulate Causal effects for SNP's based on site
            causal_eff = rng.normal(0, 1, size=(nCausal, self.nsites))

            for i in range(self.nsubjects):
                # determine the persons site
                site = self.df["abcd_site"][i]

                # multiply each persons genotype with the snp effects for the specified site
                Gene_contrib[i, ] = (
                    np.dot(Xcausal[i, :], causal_eff[:, site]))

        else:
            # sim effect from each SNP (Sample from N(0,1), later rescale to get desired variance contributions)
            causal_eff = rng.normal(0, 1, (nCausal, 1))
            Gene_contrib = np.array(Xcausal * causal_eff).flatten()

        # genetic contribution
        self.df["Gene_contrib"] = Gene_contrib

    def sim_pheno(self, var_comps=[0.5, 0.25, 0.25], phen = 1):
        # make sure no zero variance components
        for i, v in enumerate(var_comps):
            if v == 0:
                var_comps[i] = 1e-6

        # Sim errors
        errors = rng.normal(0, 1, size=self.nsubjects)

        # Calculate desired variance ratios
        StoG = var_comps[1]/var_comps[0]
        EtoG = var_comps[2]/var_comps[0]

        # Calculate empricial variance ratios
        gen_var = np.var(self.df["Gene_contrib"])
        site_var = np.var(self.df["Site_contrib"])

        # Find the empirical ratio of variances
        StoG_sim = site_var / gen_var
        site_variance_scaling = StoG_sim / StoG

        EtoG_sim = np.var(errors) / gen_var
        error_variance_scaling = EtoG_sim / EtoG

        # Scale site effects so total contributed variance is as presecribed
        self.df["Site_contrib"] = self.df["Site_contrib"] / \
            np.sqrt(site_variance_scaling)
        self.site_var = np.var(self.df["Site_contrib"])

        # Sample errors from normal scaled by the ratio of the intedended variance between genetic and error effects
        self.df["errors"] = (
            errors / np.sqrt(error_variance_scaling)).flatten()
        # Sim phenotype
        if phen==0 :
            phenoname = "Y"
        else :
            phenoname = "Y" + str(phen)
        self.df[phenoname] = self.df.Gene_contrib + self.df.Site_contrib + self.df.errors
        self.df[phenoname] = self.df[phenoname] - np.mean(self.df[phenoname])


            
    def full_sim(self, sigma, site_comp="IID",
                        nsites=30, theta_alleles=0.5, nclusts=5, dominance=5,
                        prop_causal=0.25, site_dep=False, nsubjects=1000,
                        nnpc=0, phens = 2):
        # Run through full simulation and estimation
        self.sim_sites(nsites= nsites, eq_sites=False)
        self.sim_pops(theta_alleles=0.5, nclusts=nclusts, site_comp= site_comp, dominance=2)
        self.sim_genos()
        self.sim_gen_effects(prop_causal=prop_causal, site_dep= site_dep)
        
        for i in range(phens) :
            self.sim_pheno(var_comps=sigma, phen = i)

        # Make estimates

        # # Fit basic AdjHE with RV
        # self.estimate(nnpc=nnpc, Method="AdjHE", RV="abcd_site")
        # Site_RE = self.result["h2"][0]

        # # Fit basic AdjHE
        # self.estimate(nnpc=nnpc, Method="AdjHE")
        # Basic_est = self.result["h2"][0]

        # # Fit AdjHE with Site fixed effect
        # self.estimate(nnpc=nnpc, Method="AdjHE", covars=["abcd_site"])
        # Site_FE = self.result["h2"][0]

        # # Fit AdjHE naively pooling between sites
        # self.estimate(nnpc=nnpc, Method="Naive")
        # Naive = self.result["h2"]

        # Fit MOM
        #sim1.estimate(nnpc = 0, fast = False, covars= ["abcd_site"])
        # MOM = sim1.result["h2"][0]
        
        # Fit GCTA
        # self.estimate(Method="GCTA")
        # GCTA_result = self.result["h2"]
        
        # # return results
        # self.result = pd.DataFrame({"sg" : sigma[0], "ss" : sigma[1], "se" : sigma[2],                  # Store variance parameters
        #               "nsub" : self.nsubjects, "nsites" : nsites, "nclusts" : nclusts,            # Store count parameters
        #               "prop_causal" : prop_causal, "site_comp" : site_comp, "nSNPs" : self.nSNPs, # Store other parameters
        #               "theta_alleles" : theta_alleles, "dominance" : dominance, "npc" : nnpc,
        #               "Site_RE" : Site_RE, "Basic_est" : Basic_est, "Site_FE" :Site_FE, "Naive_AdjHE" : Naive, "GCTA" : GCTA_result # Store estimates
        #               }, index = [0])
        