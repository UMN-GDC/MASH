# make sure your working directory is MASH
from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from AdjHE.estimation.all_estimators import Basu_estimation
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest

rng = np.random.default_rng(123)
# Instantiate a simulator
sim = pheno_simulator(rng = rng, nsubjects= 1000)
# Detemrine how many sites, genetic clusters you want
sim.sim_sites(nsites =1)
sim.sim_pops(nclusts= 1)
sim.sim_genos()
# Simulate a phenotype with the desired heritability coming that are homogeneous across
# clusters (h2Hom), and different between clusters (h2Het), as well as the power to weight
# the SNP effects as a function of their varainces (alpha)
sim.sim_pheno(h2Hom = 0.5, h2Het= [0], alpha = -1)


# In addition, if you're interested to use the built in estimation functions, instantiate an
# estimator object
est = Basu_estimation()
# specify the GRM and dataframe containing the phenotype
est.GRM = sim.GRM
est.df = sim.df
# Then estimate the heritability specifying the Method, fixed effects to control for
# Simulations automatically have a "Xc" covariate, and the number of pcs (npc)
result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "AdjHE", fixed_effects= ["Xc"])
