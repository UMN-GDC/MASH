# make sure your working directory is MASH
from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest


# The dataframe needs to have the subject ancestries contained in a column called "subj_ancestries"

rng = np.random.default_rng(123)
# Instantiate a simulator If the simulation is comprised of multiple ancestries, put a filepath containing the FID, IID in the same order as the genotpes, and a column with the ancestries
sim = pheno_simulator(rng = rng, 
                      plink_prefix="Estimate/examples/Genotype_files/geno", 
                      grmPrefix = "Estimate/examples/Input_files/grm", 
#                      covarFile = "Estimate/examples/Input_files/covar.txt",
                      subjAncestries= None)
# Detemrine how many sites, genetic clusters you want
sim.sim_pheno(h2Hom = 0.5, h2Het= [0, 0], alpha = 0)


#%%
est = Basu_estimation()
est.GRM = sim.GRM
est.df = sim.df
result = est.estimate(mpheno = ["Y0", "Y1"], npc = [0, 1], Method = "GCTA", fixed_effects= ["Xc"])


