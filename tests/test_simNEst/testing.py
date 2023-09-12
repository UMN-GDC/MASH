import numpy as np
from Estimate.estimators.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator

#%% Basic testing for single site and cluster
for i in range(5) :
    rng = np.random.default_rng(i)
    sim = pheno_simulator(rng = rng, nsubjects= 500)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    npcs = list(np.linspace(1,1,1).astype(int))

    sim.sim_pheno(h2Hom = 0, h2Het= np.repeat(0.5,3))
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    result = est.estimate(mpheno = ["Y0"], npc = npcs, Method = "GCTA", fixed_effects= ["Xc"], PC_effect = "fixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = npcs, Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "mixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = npcs, Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "random")
    print(result)
