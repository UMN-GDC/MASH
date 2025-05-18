from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import h2Estimation
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
    sim.sim_pheno(h2Hom = 0.6, h2Het= [0, 0, 0, 0, 0], alpha = 0)
    est = h2Estimation()
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
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "fixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "mixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [0,1,2], Method = "AdjHE", fixed_effects= ["Xc"], PC_effect = "random")
    print(result)

    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

#%%
@pytest.mark.C5S1_thorough
def test_thoroughC1S1(C5S1): 
    est= C5S1
   
    for cluster in est.df.subj_ancestries.unique() : 
        esti = h2Estimation()
        esti.df = est.df.query("subj_ancestries=="+ str(cluster))
        esti.GRM = est.GRM[esti.df.index, :][:, esti.df.index]
        esti.df = esti.df.reset_index()
        resulti = esti.estimate(mpheno = ["Y0"], npc = [0,1], Method = "GCTA", fixed_effects= ["Xc"], PC_effect = "fixed")
        print(resulti)
