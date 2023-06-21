from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from AdjHE.estimation.all_estimators import Basu_estimation
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest
#%%

@pytest.fixture
def singleSiteClust() :
    rng = np.random.default_rng(12345)
    sim = pheno_simulator(nsubjects= 1000)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 1)
    sim.sim_genos()
    sim.sim_gen_effects()
    sim.sim_covars(cov_effect = True)
    sim.sim_pheno(var_comps = [0.25,0.0001, 0.75])
    return sim 

@pytest.mark.simNest
def test_simNGCTA(singleSiteClust) :
    sim = singleSiteClust
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    result = est.estimate(mpheno = ["Y0"], npc = [1], Method = "GCTA", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.25) 

@pytest.mark.simNest
def test_simAdjHE(singleSiteClust):
    sim = singleSiteClust
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "AdjHE", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.25) 
