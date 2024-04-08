#! /usr/bin/env python3

"""
AdjHE estimator
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-06-06
"""

##############################################################
# Simulation runner for AdjHE package 
##############################################################
import os
import numpy as np
from Simulate.parser import get_args, read_flags

try : 
  os.chdir("/home/christian/Research/Stat_gen/tools/MASH")
  from Simulate.simulation_helpers.Sim_generator import pheno_simulator
  #from Simulate.summarizers.genotype_viz import plotClusters
  # from Simulate.summarizers.genotype_viz import plotClusters
  #os.chdir("/home/christian/Research/Integrative/prsPPMx")
except NameError :
  print("Already loaded appropriate modules")



def simulate_data(nSNPs=1000, nsubjects=500, nclusts=1, nphenos=2, shared=0.5, prop_causal=[0.25, 0.25], theta_alleles=[0.95, 0.25], h2Hom=0.8, h2Het=[0.1, 0.1], riskGroups = False):
   sim = pheno_simulator(nsubjects=nsubjects, nSNPs=nSNPs)
   sim.sim_sites()
   sim.sim_pops(nclusts=nclusts, theta_alleles=theta_alleles, shared=shared)
   sim.sim_genos()
   sim.sim_pheno(h2Hom=h2Hom, h2Het=h2Het, nphenos=nphenos, prop_causal=prop_causal, alpha=-1)
   if riskGroups :
     sim.df.Y0 = sim.df.Y0  + sim.df.riskGroups * 2
   sim.save_plink()
   
   

def main():

    #%% Get command line arguments
    # Get CL arguments and convert them to usable Python objects in a dictionary
    #args= read_flags({"argfile" : "Simulations/Sim_params/Adding_sites_clusts.json"})
    args = read_flags(get_args())
    print("Simulating with the following parameters:")
    print(args)
    # make an argparser for the simulate_data function
    print(args)
    simulate_data(nSNPs=args["nSNPs"], nsubjects=args["nsubjects"], nclusts=args["nclusts"], nphenos=args["nphenos"], shared=args["shared"], prop_causal=args["prop_causal"],
                  theta_alleles=args["theta_alleles"], h2Hom=args["h2Hom"], h2Het=args["h2Het"], riskGroups = args["riskGroups"])
    
    
    # sim_n_est(nsubjects = args["nsubjects"], h2 = 0.5, nsites = 30,
    #               nclusts =1, nnpc = 1,
    #               nSNPs=20, phens = 2,
    #               random_BS=True, savePlink = args["savePlink"],
    #               estimateHeritability = args["estimateHeritability"])

    
    
