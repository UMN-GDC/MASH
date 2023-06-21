from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from AdjHE.estimation.all_estimators import Basu_estimation
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest
#%%

@pytest.fixture
def singleSiteClust() :
    sim = pheno_simulator(nsubjects= 250)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 1)
    sim.sim_genos()
    sim.sim_gen_effects()
    sim.sim_covars(cov_effect = True)
    sim.sim_pheno(var_comps = [0.5,0.0001, 0.5])
    return sim 

@pytest.mark.simNGCTA
def test_simNGCTA(singleSiteClust) :
    sim = singleSiteClust
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    est.estimate(mpheno = ["Y0"], npc = [0], Method = "GCTA", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.4) 

@pytest.mark.simNAdjHE
def test_simAdjHE(singleSiteClust):
    sim = singleSiteClust
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    est.estimate(mpheno = ["Y0"], npc = [0], Method = "AdjHE", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.4) 
