import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator

reps = 100
h2s = pd.DataFrame(np.zeros((2 * 2* reps, 2)), columns = ["Combat", "SWD"]) 
h2s.index = pd.MultiIndex.from_product([["S2", "S25"], ["EQUAL", "UNEQUAL"], list(range(reps))], names= ["Sites", "Distribution", "Rep"]) 
rng = np.random.default_rng(1234)
for sites, dist, rep in itertools.product([2, 25], ["EQUAL", "UNEQUAL"], list(range(reps))) : 
    sim = pheno_simulator(rng = rng, nsubjects= 2000)
    sim.sim_sites(nsites =sites, siteDistribution = dist)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.5, h2Het= [0,0], alpha = -1)
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    h2s.loc[("S"+ str(sites), dist, rep), "Combat"] = est.estimate(mpheno = ["Y0", "Y1"], npc = [1],
                                                                 Method = "Combat", random_groups = "abcd_site", PC_effect = "random", fixed_effects= ["Xc"])["h2"][0]
    h2s.loc[("S"+ str(sites), dist, rep), "SWD"] = est.estimate(mpheno = ["Y0", "Y1"], npc = [1],
                                                                Method = "SWD", random_groups = "abcd_site", PC_effect = "random", fixed_effects= ["Xc"])["h2"][0]

h2s.to_csv("Simulate/AdjHE_sims/2k_SWD_COMBAT_0.csv", index= True)
df = pd.read_csv("Simulate/AdjHE_sims/2k_SWD_COMBAT_0.csv")
# Density plot of X
sns.boxplot(df)
# add horizontal line at 0.5
plt.axhline(0.5, color='r', linestyle='-')
plt.show()
