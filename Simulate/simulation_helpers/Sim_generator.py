#im_/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
import logging
import os
import subprocess
import numpy as np
import pandas as pd
from pandas_plink import read_plink, write_plink1_bin
from xarray import DataArray
from Estimate.data_input.load_data import ReadGRMBin
from Simulate.simulation_helpers.sites import sim_sites
from Simulate.simulation_helpers.clusters import sim_pop_alleles, assign_clusters
from Simulate.simulation_helpers.genos import sim_genos
from Simulate.simulation_helpers.pheno import sim_pheno
from Simulate.simulation_helpers.plink_pheno import sim_plink_pheno
from Estimate.estimators.gwaPRS import fitPlink2GWAS


class pheno_simulator():
    def __init__(self, rng= None, nsubjects=1000, nSNPs=1000, plink_prefix=None, grmPrefix= None, covarFile=  None, subjAncestries = None):
        self.rng = np.random.default_rng(rng)
        self.plink_prefix = plink_prefix
        
        if plink_prefix is None :
            self.nsubjects = nsubjects
            self.nSNPs = nSNPs

            # rng the simulated dataframe
            self.df = pd.DataFrame({"FID": np.arange(start=1, stop=nsubjects + 1), 
                                    "IID": np.arange(start=1, stop=nsubjects + 1)})


        else :
            if grmPrefix is None :
                logging.info("GRM prefix not specified, but that's OK")
            else : 
                self.GRM = ReadGRMBin(grmPrefix)


            self.plink_prefix= plink_prefix
            
            (bim, fam, bed) = read_plink(plink_prefix)
            (self.nSNPs, self.nsubjects) = bed.shape
            self.df = pd.DataFrame({"FID": fam.fid,
                                    "IID": fam.iid})

            # Don't need to specify sharedIdx if genotypes prespecified
            self.sharedIdx =[0] 

            if covarFile is not None :
                covars = pd.read_csv(covarFile, sep='\s+')
                self.df = pd.merge(self.df, covars, on = ["FID", "IID"])
            # bed is (nSNPs x nsubjects) so transpose it
            self.genotypes = bed.T
            self.rng = np.random.default_rng()
            # Include subject ancestries
            if subjAncestries is None :
                self.df["subj_ancestries"] = 1  
            else : 
                temp = pd.read_csv(subjAncestries, sep='\s+')
                self.df["subj_ancestries"]  = temp.iloc[:,2]


    def sim_sites(self, nsites=1, siteDistribution = "EQUAL", random_BS = True):
        self.nsites = nsites
        self.siteDistribution = siteDistribution
        self.df[["abcd_site", "Site_contrib"]] = sim_sites(rng = self.rng, nsubjects = self.nsubjects, 
                                                           nsites=nsites, siteDistribution = self.siteDistribution, 
                                                           random_BS = random_BS)

    def sim_pops(self, theta_alleles = 0.5, nclusts=1, shared = 0.5):
        self.theta_alleles = theta_alleles
        # Clusts are distinct ancestries
        self.nclusts = nclusts
        self.shared = shared
        

        # Simualte allele frequencies for common ancestor and for genetic clusters
        self.ancest_freqs, self.cluster_frequencies, self.sharedIdx = sim_pop_alleles(
            rng = self.rng,
            theta_alleles = self.theta_alleles,
            nclusts=self.nclusts, nSNPs = self.nSNPs,
            shared = shared)

        # Sample ancesrties for each individual Dim = (nsubjects,)
        self.df["subj_ancestries"] = assign_clusters(
            df = self.df,
            rng = self.rng,
            nclusts=self.nclusts,
            nsites = self.nsites)

    def sim_genos(self, maf = 0.01):
        self.maf = maf
        (self.genotypes, self.GRM, pcs) = sim_genos(rng = self.rng,
                                                    cluster_frequencies = self.cluster_frequencies, 
                                                    subject_ancestries = self.df["subj_ancestries"], maf = maf) 
        # update the number of SNPs
        self.nSNPS_post_filter = self.genotypes.shape[1]
        self.df = pd.concat([self.df, pcs], axis=1)

    def sim_pheno(self, h2Hom=0.5, h2Het=[0], alpha = -1, phenobasename = "Y", nphenos = 2, prop_causal = [0.1, 0.1]):
        self.h2Hom = h2Hom
        self.h2Het = h2Het
        self.alpha = alpha
        self.prop_causal = prop_causal
        (self.df, self.causals, self.homo_eff, self.het_eff)  = sim_pheno(rng = self.rng,
                                                                          genotypes = self.genotypes,
                                                                          df = self.df,
                                                                          h2Hom = self.h2Hom,
                                                                          h2Het = self.h2Het,
                                                                          alpha = self.alpha,
                                                                          prop_causal = prop_causal,
                                                                          sharedIdx = self.sharedIdx,
                                                                          phenoname = phenobasename + "0")
        for i in range(1, nphenos): 
            (temp, __1, __2, __3)  = sim_pheno(rng = self.rng,
                                               genotypes = self.genotypes,
                                               df = self.df,
                                               h2Hom = self.h2Hom,
                                               h2Het = self.h2Het,
                                               alpha = self.alpha,
                                               prop_causal = [0.1, 0.1],
                                               sharedIdx = self.sharedIdx,
                                               phenoname = phenobasename + str(i))
            self.df[phenobasename + str(i)] = temp[phenobasename + str(i)]
            
        self.df["riskGroups"] = np.repeat([0,1], self.df.shape[0]/ 2)
    
    def sim_outcome(self, beta):

    def full_sim(self, h2 = 0.5,
                 nsites=30, nclusts=5,
                 nsubjects=1000,
                 nnpc=0, phens = 2,
                 alpha=-1,
                 random_BS = True, npcs = 0,  maf_filter= 0.1):
        
        # Only simulate genes if plink file not specified
        if self.plink_prefix == None: 
  
            # Run through full simulation and estimation
            self.sim_sites(nsites= nsites, random_BS = random_BS)
            
            self.sim_pops(nclusts=nclusts)
            self.sim_covars()
            
            for i in range(phens) :
                self.sim_pheno(h2 = h2, phenoname = "Y" + str(i), alpha = alpha)
        
        else : 
            for i in range(phens) :
                phenoname = "Y" + str(i)
                self.df[phenoname] = sim_plink_pheno(rng = self.rng, bed = self.genotypes, h2= h2, npcs=npcs)

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

    def save_plink(self, prefix = "temp/simulation"):
        """
        Save the simulated genotypes to a plink file. By default this is a temp folder in MASH

        Parameters
        ----------
        prefix : TYPE, optional
            DESCRIPTION. The default is "sim_genostest".

        Returns
        -------
        None.

        """
        # make directory if doesn't exist
        os.makedirs(prefix, exist_ok = True)  
        
        self.prefix = prefix
        G = DataArray(
            # n x nSNPs array of genotypes
            self.genotypes.astype("float32"),
            dims=["sample", "variant"],
            coords = dict(
                sample  = self.df.IID, # IID
                fid     = ("sample", self.df.FID), #FID
                variant = np.arange(1, self.genotypes.shape[1] + 1),
                snp     = ("variant", range(1, self.genotypes.shape[1] + 1)),
                chrom   = ("variant", np.repeat(1, self.genotypes.shape[1])),
                a0      = ("variant", np.repeat("A", self.genotypes.shape[1])),
                a1      = ("variant", np.repeat("T", self.genotypes.shape[1])),
            )
        
        )

        write_plink1_bin(G, f"{prefix}.bed")

        # Save phenotype separate from subj_ancestries and covariates
        self.df.to_csv(f"{prefix}.covar", sep = "\t", index = False)
       
        # Save FID, ID and columns starting with Y
        self.df["P1"] = 0
        self.df["P2"] = 0
        self.df["sex"] = 0
        self.df[["FID", "IID", "P1", "P2", "sex", "sex"]].to_csv(f"{prefix}.fam", sep = "\t", index = False, header=False)
        self.df[["FID", "IID"] + list(self.df.filter(regex='^Y'))].to_csv(f"{prefix}.pheno", index = False, header=True, sep = "\t")

        # Save subject_ancestries
        try : 
            self.df[["FID", "IID", "subj_ancestries", "abcd_site", "riskGroups", "confound"]].to_csv(f"{prefix}.covar", sep = "\t", index = False)
        except :
            self.df[["FID", "IID", "subj_ancestries", "abcd_site", "riskGroups"]].to_csv(f"{prefix}.covar", sep = "\t", index = False)

        # Read the output .bim file
        bim = pd.read_table(f"{prefix}.bim", sep='\s+', header = None)
        # replace the thrid and fourth columns with integers going from 1 to nSNPs
        bim.iloc[:,2] = np.arange(self.genotypes.shape[1])

        bim.iloc[:,3] = np.arange(self.genotypes.shape[1])
        # rewrite the table 
        bim.to_csv(f"{prefix}.bim", sep = "\t", index = False, header = False)
        
       # command = f"plink --bfile {prefix} --r2 --no-parents --out {prefix}"
       # process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
       # 
       # # Fit PRS
       # command = f"plink2 --bfile {prefix} --glm --covar {prefix}.covar --out {prefix}"
       # process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
        
        
        
    