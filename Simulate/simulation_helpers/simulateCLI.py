import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from xarray import DataArray
import pandas_plink as pdplink
import os
import argparse

try : 
  os.chdir("/home/christian/Research/Stat_gen/tools/MASH")
  from Simulate.simulation_helpers.Sim_generator import pheno_simulator
  from Simulate.summarizers.genotype_viz import plotClusters
  # from Simulate.summarizers.genotype_viz import plotClusters
  os.chdir("/home/christian/Research/Integrative/prsPPMx")
except NameError :
  print("Already loaded appropriate modules")
  
  
def simulate_data(nSNPs=1000, nsubjects=500, nclusts=1, nphenos=2, shared=0.5, prop_causal=[0.25, 0.25], theta_alleles=[0.95, 0.25], h2Hom=0.8, h2Het=[0.1, 0.1]):
   sim = pheno_simulator(nsubjects=nsubjects, nSNPs=nSNPs)
   sim.sim_sites()
   sim.sim_pops(nclusts=nclusts, theta_alleles=theta_alleles, shared=shared)
   sim.sim_genos()
   sim.sim_pheno(h2Hom=h2Hom, h2Het=h2Het, nphenos=nphenos, prop_causal=prop_causal, alpha=-1)
   sim.save_plink()
   
   
# make an argparser for the simulate_data function
#print(args)
#simulate_data(nSNPs=args.nSNPs, nsubjects=args.nsubjects, nclusts=args.nclusts, nphenos=args.nphenos, shared=args.shared, prop_causal=args.prop_causal, theta_alleles=args.theta_alleles, h2Hom=args.h2Hom, h2Het=args.h2Het)
