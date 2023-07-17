from Estimate.estimation.GCTA_wrapper import GCTA, gcta
from Estimate.estimation.all_estimators import Basu_estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest

#%% Basic testing for single site and cluster

@pytest.fixture
def C5S1() :
    rng = np.random.default_rng(123)
    sim = pheno_simulator(rng = rng, nsubjects= 1000)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 5)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.2, h2Het= [0.2, 0.4, 0.6, 0.8, 0.0], alpha = -1)
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    return est 

@pytest.mark.C5S1
def test_simNGCTA(C5S1) :
    est = C5S1
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "GCTA", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

@pytest.mark.C5S1
def test_simAdjHE(C5S1):
    est = C5S1
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "AdjHE", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

#%%
@pytest.mark.C5S1_thorough
def test_thoroughC1S1(C5S1): 
    est= C5S1
   
    for cluster in est.df.subj_ancestries.unique() : 
        esti = Basu_estimation()
        esti.df = est.df.query("subj_ancestries=="+ str(cluster))
        esti.GRM = est.GRM[esti.df.index, :][:, esti.df.index]
        resulti = esti.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "GCTA", fixed_effects= ["Xc"])
        print(resulti)
