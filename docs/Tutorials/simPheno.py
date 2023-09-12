from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np

rng = np.random.default_rng(123)
sim = pheno_simulator(rng = rng, nsubjects= 1000)
sim.sim_sites(nsites =1)
sim.sim_pops(nclusts= 2)
sim.sim_genos()
sim.sim_pheno(h2Hom = 0.5, h2Het= [0, 0], alpha = 0)

sim.df.to_csv('simulated_pheno.csv', index=False)

