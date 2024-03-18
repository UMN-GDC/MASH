from ctypes import CDLL 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from xarray import DataArray
import pandas_plink as pdplink
import os
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
from Simulate.summarizers.genotype_viz import plotClusters

nSNPs = 1000
nsubjects = 500
nclusts = 1
nphenos = 2
shared = 0.5
prop_causal = [0.25, 0.25]
theta_alleles = [0.95, 0.25]
h2Hom =0.8
# Add this for transfer learning 
h2Het = [0.1, 0.1]

rng = np.random.default_rng(12)
nclusts = 2 


sim = pheno_simulator(nsubjects = nsubjects, nSNPs = nSNPs)
sim.sim_sites()
sim.sim_pops(nclusts = nclusts, theta_alleles = theta_alleles, shared = shared)
sim.sim_genos()
sim.sim_pheno(h2Hom =h2Hom, h2Het = h2Het, nphenos = nphenos, prop_causal = prop_causal, alpha =-1)
# Subtract the mean from every column starting with Y and store them in the dataframe
sim.df.iloc[:,30:(30+nphenos)] = sim.df.filter(regex='^Y') - sim.df.filter(regex= "^Y").mean()

sim.save_plink("test2")
sim.fitGWAS(prefix = "test2")


