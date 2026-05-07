import itertools
import numpy as np
import pandas as pd
from Estimate.estimators.all_estimators import h2Estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator

reps = 50
h2s = pd.DataFrame(np.zeros((2 * 2* reps, 2)), columns = ["Combat", "SWD"]) 
h2s.index = pd.MultiIndex.from_product([["S2", "S25"], ["EQUAL", "UNEQUAL"], list(range(reps))], names= ["Sites", "Distribution", "Rep"]) 
rng = np.random.default_rng(1234)
for sites, dist, rep in itertools.product([2, 25], ["EQUAL", "UNEQUAL"], list(range(reps))) : 
    sim = pheno_simulator(rng = rng, nsubjects= 5000)
    sim.sim_sites(nsites =sites, siteDistribution = dist)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.66, h2Het= [0,0], alpha = -1)
    est = h2Estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    # Set args for the estimate method
    est.args = {
        "qcovar": ["Xc"],
        "covar_discrete": [],
        "npc": [1],
        "Method": None,  # Will be set before each call
        "loop_covars": False,
        "Naive": False,
        "random_groups": "abcd_site",
        "covar": None,
    }
    est.args["Method"] = "Combat"
    h2s.loc[("S"+ str(sites), dist, rep), "Combat"] = est.estimate(PC_effect = "random")["h2"][0]
    est.args["Method"] = "SWD"
    h2s.loc[("S"+ str(sites), dist, rep), "SWD"] = est.estimate(PC_effect = "random")["h2"][0]

h2s.to_csv("Simulate/AdjHE_sims/5k_SWD_COMBAT_0.csv", index= True)

reps = 50
h2s = pd.DataFrame(np.zeros((2 * 2* reps, 2)), columns = ["Combat", "SWD"]) 
h2s.index = pd.MultiIndex.from_product([["S2", "S25"], ["EQUAL", "UNEQUAL"], list(range(reps))], names= ["Sites", "Distribution", "Rep"]) 
rng = np.random.default_rng(1234)
for sites, dist, rep in itertools.product([2, 25], ["EQUAL", "UNEQUAL"], list(range(reps))) : 
    sim = pheno_simulator(rng = rng, nsubjects= 5000)
    sim.sim_sites(nsites =sites, siteDistribution = dist)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.66, h2Het= [0,0], alpha = -1)
    est = h2Estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    est.args = {
        "qcovar": ["Xc"],
        "covar_discrete": [],
        "npc": [1],
        "Method": None,  # Will be set before each call
        "loop_covars": False,
        "Naive": False,
        "random_groups": "abcd_site",
        "covar": None,
    }
    est.args["Method"] = "Combat"
    h2s.loc[("S"+ str(sites), dist, rep), "Combat"] = est.estimate(PC_effect = "random")["h2"][0]
    est.args["Method"] = "SWD"
    h2s.loc[("S"+ str(sites), dist, rep), "SWD"] = est.estimate(PC_effect = "random")["h2"][0]
h2s.to_csv("Simulate/AdjHE_sims/5k_SWD_COMBAT_1.csv", index= True)

