from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import h2Estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest

#%% Basic testing for single site and cluster
for i in range(10) : 
    rng = np.random.default_rng(i)
    sim = pheno_simulator(rng = rng, nsubjects= 1000)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 5)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.6, h2Het= [0, 0, 0, 0, 0], alpha = 0)
    est = h2Estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    result = est.estimate(mpheno = ["Y0"], npc = [1,2,3], Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "fixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [1,2,3], Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "mixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [1,2,3], Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "random")
    print(result)
