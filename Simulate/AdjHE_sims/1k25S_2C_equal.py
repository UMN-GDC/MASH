from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
reps = 100 
h2s = np.zeros(reps)
rng = np.random.default_rng(123)


for i in range(reps) : 
    sim = pheno_simulator(rng = rng, nsubjects= 1000)
    sim.sim_sites(nsites =25, siteDistribution = "EQUAL")
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.5, h2Het= [0,0], alpha = -1)
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    h2s[i]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "COMBAT", random_groups = "abcd_site", PC_effect = "random", fixed_effects= ["Xc"])["h2"][0]


np.save("Simulate/AdjHE_sims/1k25S_2C_equal.npy", h2s)

#%%
import matplotlib.pyplot as plt
import seaborn as sns

X = np.load("Simulate/AdjHE_sims/1k25S_2C_equal.npy")

# Density plot of X
sns.boxplot(X)
# add horizontal line at 0.5
plt.axhline(0.5, color='r', linestyle='-')
plt.show()
