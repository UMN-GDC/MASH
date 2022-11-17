#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
from plotly.offline import plot
import plotly.express as px
import itertools
import os
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")

import numpy as np
import pandas as pd
from functions.Estimation.all_estimators import load_n_estimate
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
rng = np.random.default_rng(123)


# TO DO:
# Compare with naively splitting into sites and estimating (see if heritability estimate changes by site)


# %%
class AdjHE_simulator():
    def __init__(self, nsubjects, nSNPs):
        self.nsubjects = nsubjects
        self.nSNPs = nSNPs

        # Seed the dataframe
        self.df = pd.DataFrame({"fid": np.arange(start=1, stop=nsubjects + 1),
                                "iid": np.arange(start=1, stop=nsubjects + 1),
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

    def pop_clusts(self, npc=2):
        # Checked genotypes with coming from separate clusters
        trans = PCA(n_components=npc).fit_transform(self.genotypes)
        sns.scatterplot(x=trans[:, 0], y=trans[:, 1],
                        hue=self.df.subj_ancestries)

    def sim_gen_effects(self, prop_causal=0.05, site_dep=False):
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

    def sim_pheno(self, var_comps=[0.5, 0.25, 0.25]):
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
        self.df["Y"] = self.df.Gene_contrib + \
            self.df.Site_contrib + self.df.errors
        self.df["Y"] = self.df["Y"] - np.mean(self.df["Y"])

    def GRM_vis(self, sort_by=None, location=None, npc=0):
        if sort_by == None:
            df = self.df
        else:
            df = self.df.sort_values(sort_by)

        # Arrange the GRM in the same manner
        G = self.GRM[df.index, :][:, df.index]

        if npc != 0:
            # subtract substructre from GRM
            pcs = np.matrix(PCA(n_components=npc).fit(self.GRM).components_.T)
            P = pcs * np.linalg.inv(pcs.T * pcs) * pcs.T
            G = (np.eye(self.nsubjects) - P) * self.GRM

        plt.imshow(G, )
        # major_ticks = np.unique(df.abcd_site, return_counts= True)[1].cumsum()
        # plt.xticks(major_ticks -1)
        # plt.yticks(np.flip(major_ticks))
        # plt.grid(color='k', linestyle='-', linewidth=0.5)
        plt.axis('off')

        plt.title(("GRM: clusters = {c}, Site pops = {site}"
                   .format(c=self.nclusts, site=self.site_comp)))
        plt.colorbar()

        if location != None:
            plt.savefig('docs/Presentations/Images/GRM_{n}_{s}_{c}_{site}.png'.format(
                n=self.nsubjects, s=self.nsites, c=self.nclusts, site=self.site_comp),
                dpi=200)

    def estimate(self, Method=None, RV=None, nnpc=0, covars=[], silent=True):

        if (nnpc != 0) and (RV != None):
            # project data onto pc space
            pcs = pd.DataFrame(PCA(n_components=nnpc).fit_transform(self.GRM))
            pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
            # add the pcs to the dataframe
            self.df = pd.concat([self.df, pcs], axis=1)

            # subtract substructre from GRM
            pcs = np.matrix(PCA(n_components=nnpc).fit(self.GRM).components_.T)
            P = pcs * np.linalg.inv(pcs.T * pcs) * pcs.T
            G = (np.eye(self.nsubjects) - P) * self.GRM
        else:
            G = self.GRM

        if Method == "GCTA":
            # Get lower triangle indices
            l = np.tril_indices(self.nsubjects)

            # Save the GRM and phenotype for use
            self.GRM[l].astype('f4').tofile("temp.grm.bin")
            self.df[["fid", "iid", "Y"]].to_csv(
                "temp.phen", sep=" ", index=False, header=False)
            self.df[["fid", "iid"]].to_csv(
                "temp.grm.id", sep=" ", index=False, header=False)

            # Estimate using GCTA
            bashcommand = "../GCTA/gcta_1.93.2beta/gcta64 --grm temp --pheno temp.phen --mpheno 1 --reml --out temp"
            process = subprocess.Popen(
                bashcommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()

            # read the GCTA results
            df = pd.read_table("temp.hsq", sep="\t")

            self.result = {"h2": df["Variance"][df.Source == "V(G)/Vp"].item(),
                           "SE": df["SE"][df.Source == "V(G)/Vp"].item(),
                           "Pheno": "Y",
                           "PCs": nnpc,
                           "Covariates": "NONE",
                           "Time for analysis(s)": 0,
                           "Memory Usage": 0}

        elif Method == "AdjHE":
            self.result = load_n_estimate(df=self.df, covars=covars,  nnpc=nnpc, mp="Y", GRM=G,
                                          std=False, Method=Method, RV=RV, silent=silent)

        elif Method == "Naive":
            # Empty results list
            results = pd.DataFrame({"Estimate": [],
                                   "Size": []})

            # loop over  all sites
            for site in np.unique(self.df.abcd_site):
                # Grab the portion that lies within a given site
                subdf = self.df.loc[self.df.abcd_site == site, :]
                # Get size
                size = subdf.shape[0]
                # Estimate just on the supsample
                result = load_n_estimate(df=subdf, covars=[],  nnpc=0, mp="Y", GRM=self.GRM, std=False, Method="AdjHE", RV=None,
                                         silent=silent)
                result = pd.DataFrame({"Estimate": [result["h2"][0]],
                                       "Size": [size]})
                # Add to the list of estimates
                results = results.append(result, ignore_index=True)

            # Pool the estimates
            results["nXbar"] = (
                results["Size"] * results["Estimate"]) / self.nsubjects
            final_result = np.sum(results["nXbar"])
            self.result = {"h2": final_result,
                           "SE": 0,
                           "Pheno": "Y",
                           "PCs": nnpc,
                           "Covariates": "NONE",
                           "Time for analysis(s)": 0,
                           "Memory Usage": 0}
            
    def sim_n_est(self, sigma, site_comp="IID",
                        nsites=30, theta_alleles=0.5, nclusts=5, dominance=5,
                        prop_causal=0.7, site_dep=False,
                        nnpc=0):
        # Run through full simulation and estimation
        self.sim_sites(nsites= nsites, eq_sites=False)
        self.sim_pops(theta_alleles=0.5, nclusts=nclusts, site_comp= site_comp, dominance=2)
        self.sim_genos()
        self.sim_gen_effects(prop_causal=prop_causal, site_dep= site_dep)
        self.sim_pheno(var_comps=sigma)

        # Make estimates

        # Fit basic AdjHE with RV
        self.estimate(nnpc=nnpc, Method="AdjHE", RV="abcd_site")
        Site_RE = self.result["h2"][0]

        # Fit basic AdjHE
        self.estimate(nnpc=nnpc, Method="AdjHE")
        Basic_est = self.result["h2"][0]

        # Fit AdjHE with Site fixed effect
        self.estimate(nnpc=nnpc, Method="AdjHE", covars=["abcd_site"])
        Site_FE = self.result["h2"][0]

        # Fit AdjHE naively pooling between sites
        self.estimate(nnpc=nnpc, Method="Naive")
        Naive = self.result["h2"]

        # Fit MOM
        #sim1.estimate(nnpc = 0, fast = False, covars= ["abcd_site"])
        # MOM = sim1.result["h2"][0]
        
        # Fit GCTA
        self.estimate(Method="GCTA")
        GCTA = self.result["h2"]
        
        # return results
        self.result = pd.DataFrame({"sg" : sigma[0], "ss" : sigma[1], "se" : sigma[2],                  # Store variance parameters
                      "nsub" : self.nsubjects, "nsites" : nsites, "nclusts" : nclusts,            # Store count parameters
                      "prop_causal" : prop_causal, "site_comp" : site_comp, "nSNPs" : self.nSNPs, # Store other parameters
                      "theta_alleles" : theta_alleles, dominance : "dominance", "npc" : nnpc,
                      "Site_RE" : Site_RE, "Basic_est" : Basic_est, "Site_FE" :Site_FE, "Naive_AdjHE" : Naive, "GCTA" : GCTA # Store estimates
                      }, index = [0])
        
def sim_experiment(nsubjectss = [100], site_comps=["IID"], nSNPss = [20],
                    nsitess=[30], theta_alleless=[0.5], nclustss=[5], dominances=[5],
                    prop_causals=[0.7], site_deps=[False], reps=25,
                    nnpcs=[0], sigmas = [[0.5,0.25,0.25]]) :
    # Seed empty dataframe
    sim_results = pd.DataFrame()
    # loop over all simulation variables and number of repetitions
    for nsubjects, sigma, site_comp, nsites, theta_alleles, nclusts, dominance, prop_causal, site_dep, nnpc, __, nSNPs in itertools.product(nsubjectss, sigmas, site_comps, nsitess,
                                                                                                                theta_alleless, nclustss, dominances, 
                                                                                                                prop_causals, site_deps, 
                                                                                                                nnpcs, range(reps), nSNPss) :
        sim = AdjHE_simulator(nsubjects = nsubjects, nSNPs = nSNPs)
        sim.sim_n_est(sigma = sigma, site_comp = site_comp,
                            nsites= nsites, theta_alleles = theta_alleles, nclusts = nclusts, dominance= dominance,
                            prop_causal=prop_causal, site_dep=site_dep, nnpc=0)
        
        sim_results = pd.concat([sim_results, sim.result], ignore_index = True)
        
    return sim_results




def plot_save(sim_results, out) :

