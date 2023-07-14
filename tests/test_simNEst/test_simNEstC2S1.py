from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from AdjHE.estimation.all_estimators import Basu_estimation
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest

#%% Basic testing for single site and cluster

@pytest.fixture
def C2S1() :
    rng = np.random.default_rng(123)
    sim = pheno_simulator(rng = rng, nsubjects= 500)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.4, h2Het= [0.1, 0.1], alpha = -1)
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    return est 

@pytest.mark.C2S1
def test_simNGCTA(C2S1) :
    est = C2S1
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "GCTA", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

@pytest.mark.C2S1
def test_simAdjHE(C2S1):
    est = C2S1
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "AdjHE", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

#%%
@pytest.mark.C1S1_thorough
def test_thoroughC1S1(): 
    rng = np.random.default_rng(12)
    h2s = np.zeros(25)
    for i in range(25) : 
        sim = pheno_simulator(rng = rng, nsubjects= 250)
        sim.sim_sites(nsites =1)
        sim.sim_pops(nclusts= 1)
        sim.sim_genos()
        sim.sim_gen_effects()
        sim.sim_covars()
        sim.sim_pheno(h2 = 0.5)
        est = Basu_estimation()
        est.GRM = sim.GRM
        est.df = sim.df
        h2s[i]= est.estimate(mpheno = ["Y0"], npc = [0], Method = "GCTA", fixed_effects= ["Xc"])["h2"][0]
