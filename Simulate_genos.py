#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
import os 
#os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")

import numpy as np
import pandas as pd
from functions.AdjHE_estimator import load_n_estimate
from sklearn.decomposition import PCA
rng = np.random.default_rng(12345)


#%%
class AdjHE_simulator() :
    def __init__(self, nsubjects, nSNPs):
        self.nsubjects= nsubjects
        self.nSNPs = nSNPs

               
    def sim_pops(self, theta_alleles = 0.5, nclusts =1) :
        # Theta alleles controls how much ancstry is conserved. low theta_alleles = high conservation
        self.theta_alleles = theta_alleles
        # Clusts are distinct ancestries
        self.nclusts = nclusts
        
        # simulate ancestral frequency of each SNP, dimensions = (SNPs,)
        self.ancest_freqs = rng.uniform(low = 0.1, high = 0.9, size = self.nSNPs)
        
        #simulate allele frequency for each population, Dim = (nclusts x SNPS)
        self.pop_freqs = rng.beta(self.ancest_freqs * (1-self.theta_alleles)/ self.theta_alleles,
                                      (1-self.ancest_freqs) * (1-self.theta_alleles)/self.theta_alleles, size = (self.nclusts, self.nSNPs))
        
        # Sample ancesrties for each individual Dim = (nsubjects,)
        self.subj_ancestries = rng.integers(low = 0, high = self.nclusts, size = self.nsubjects)
        
        # Checked samples came from ancester with the following
        # plt.scatter(sim1.ancest_freqs, sim1.pop_freqs.mean(axis = 0))

        
    def sim_genos(self) :
        # simulate genotypes
        genotypes = rng.binomial(n = np.repeat(2, self.nSNPs), p = self.pop_freqs[self.subj_ancestries], size = (self.nsubjects, self.nSNPs))
        # keep SNPs with MAF greater than 0.05
        maf_filter = np.logical_and((np.sum(genotypes, axis = 0) / (2* self.nsubjects))> 0.05,
                                    (np.sum(genotypes, axis = 0) / (2* self.nsubjects))< 0.95)
        self.genotypes = genotypes[:,  maf_filter]
        
        
        # get new number of SNPs
        self.nSNPs = self.genotypes.shape[1]
        # Checked genotypes with coming from separate clusters
        # trans = PCA(n_components=2).fit_transform(sim1.genotypes)
        # sns.scatterplot(x= trans[:,0], y= trans[:,1], hue = sim1.subj_ancestries)

        # reference allele frequ
        # standardize the genotpyes 
        allele_freqs = np.mean(self.genotypes, axis = 0) /2 
        genotypes = np.matrix((self.genotypes - 2 * allele_freqs) / np.sqrt(2 * allele_freqs * (1- allele_freqs)))
        self.GRM = np.dot(genotypes, genotypes.T) / self.nSNPs



    def sim_phenos(self, var_comps = [0.5,0, 0.5], intercept = 0, nsites =1, prop_causal=0.05) :
        # make sure no zero variance components
        for i, v in enumerate(var_comps) :
            if v == 0:
                var_comps[i] = 1e-6
        # number of locations
        self.nsites = nsites
        self.prop_causal = prop_causal
        nSNPs = self.nSNPs        
        nCausal =  int(nSNPs * prop_causal)
        # select causal snos
        causals = rng.choice(nSNPs, nCausal, replace = False, shuffle = False)
        
        # sim effect from each SNP
        causal_eff = rng.normal(np.repeat(0, nCausal), np.repeat(var_comps[0] /nCausal, nCausal))

        # seed empty effect vector
        all_eff = np.zeros(nSNPs)
        # input simulated causal SNP's
        all_eff[causals] = causal_eff
        # return the effects vector
        all_eff = np.matrix(all_eff).T

        # genetic contribution
        Xu = self.genotypes * all_eff
        gen_var = np.var(Xu)
        self.gen_var = gen_var
        
        # Calculate Site variance
        Groups = np.repeat(np.arange(start = 1, stop = nsites + 1), int(self.nsubjects/nsites))
        # Make a dummy matrix
        Groups_dum = np.matrix(pd.get_dummies(Groups))
        ## standardize it 
        # group_freq = np.mean(Groups_dum, axis = 0)
        # Groups_dum = np.matrix((Groups_dum - group_freq) / np.sqrt(2 * group_freq * (1- group_freq)))
        StoG = var_comps[1]/var_comps[0]
        Bs = np.matrix(np.random.normal(0,  1, nsites)).T
        # find n per group
        site_effects = Groups_dum * Bs
        site_var = np.var(site_effects)
        # Scale site effects so total contributed variance is as presecribed
        StoG_sim = site_var / gen_var
        site_variance_scaling = StoG_sim / StoG
        self.site_variance_scaling = site_variance_scaling
        site_effects = site_effects / np.sqrt(site_variance_scaling)      
        self.site_var = np.var(site_effects)
        
        # Sample errors from normal scaled by the ratio of the intedended variance between genetic and error effects
        errors = np.matrix(rng.normal(0, 1, size = self.nsubjects)).T
        EtoG = var_comps[2]/var_comps[0]
        EtoG_sim = np.var(errors) / gen_var
        error_variance_scaling = EtoG_sim / EtoG
        errors = errors / np.sqrt(error_variance_scaling)      
        self.err_var= np.var(errors)

        
        df = pd.DataFrame({"fid" : Groups,
                           "iid" : np.arange(start= 1, stop = self.nsubjects + 1),
                           "abcd_site" : Groups,
                           "Y" : list(np.array(Xu + errors + intercept + site_effects).flatten())})


        # return simulated Y and genotypes
        self.df = df
        
        
    def estimate(self, fast = True, RV = None, nnpc = 0, covars = []):
        
        if nnpc != 0 :
            # project data onto pc space
            pcs = pd.DataFrame(PCA(n_components=nnpc).fit_transform(self.GRM))
            pcs.columns = ["pc_" + str(col + 1) for col in pcs.columns]
            # add the pcs to the dataframe
            self.df = pd.concat([self.df, pcs], axis = 1)
        
        
        self.result = load_n_estimate(df=self.df, covars=covars,  nnpc=nnpc, mp="Y", GRM= self.GRM, std= False, fast = fast, RV = RV)


#%% Simulations
import itertools
#import plotly.express as px
#from plotly.offline import plot

sg = [0.25, 0.5, 0.75]
ss = [0, 0.5]
se = [0, 0.25, 0.5, 0.75]

sigmas = []
for ses in itertools.product(sg,ss, se) :
    if (sum(ses) == 1) and (ses[0] > 0.1) :
        sigmas.append(ses)
      
# Add a completely null case
sigmas.insert(0, (0, 0, 1))

#%%
results = {"sg" : [], 
           "ss" : [],
           "se" : [],
           "Basic_est" : [],
           "Site_RE" : [],
           "Site_FE" : []}
for sigma in sigmas :
    for rep in range(100) :
        results["sg"].append(sigma[0])
        results["ss"].append(sigma[1])
        results["se"].append(sigma[2])
    
        
        sim1 = AdjHE_simulator(nsubjects= 2000, nSNPs = 100)
        sim1.sim_pops(theta_alleles = 0.5, nclusts = 1)
        sim1.sim_genos()
        sim1.sim_phenos(prop_causal = 0.05, var_comps= list(sigma), nsites=  25, intercept = 0)
        # Fit basic AdjHE
        sim1.estimate(nnpc = 0, fast = True, RV = "abcd_site")
        results["Site_RE"].append(sim1.result["h2"][0])
    
        # Fit AdjHE with Site RV
        sim1.estimate(nnpc = 0, fast = True)
        results["Basic_est"].append(sim1.result["h2"][0])
        
        # Fit AdjHE with Site RV
        sim1.estimate(nnpc = 0, fast = True, covars= ["abcd_site"])
        results["Site_FE"].append(sim1.result["h2"][0])

        
results = pd.DataFrame(results)
results2 = pd.melt(results, id_vars= ["sg", "ss", "se"], value_vars=['Basic_est', 'Site_RE', 'Site_FE'])
results2["herit"] = results2.sg/ (results2.sg + results2.ss + results2.se)
results.to_csv("sim_2000.csv",index = False)

#%%


#fig = px.box(results2, x="variable", y="value", facet_col="herit", facet_row="ss" )
#fig.add_hline(y = "sg")
# Free y scals
#fig.update_yaxes(matches=None)
#fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))

# Save figure
#plot(fig, filename='full_genome_1000_sim_results.html')
