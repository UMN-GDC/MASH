from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from AdjHE.estimation.all_estimators import Basu_estimation
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest
#%%

@pytest.fixture
def singleSiteClust() :
    rng = np.random.default_rng(123)
    sim = pheno_simulator(rng = rng, nsubjects= 1000)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 1)
    sim.sim_genos()
    sim.sim_gen_effects()
    sim.sim_covars()
    sim.sim_pheno(h2 = 0.5)
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    return est 

@pytest.mark.simNest
def test_simNGCTA(singleSiteClust) :
    est = singleSiteClust
    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "GCTA", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

@pytest.mark.simNest
def test_simAdjHE(singleSiteClust):
    est = singleSiteClust
    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "AdjHE", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 
