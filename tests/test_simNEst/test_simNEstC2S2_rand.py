from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import h2Estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest

#%% Basic testing for single site and cluster

@pytest.fixture
def C2S2_rand() :
    rng = np.random.default_rng(123)
    sim = pheno_simulator(rng = rng, nsubjects= 1000)
    sim.sim_sites(nsites =2, siteDistribution = "EQUAL")
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    sim.sim_pheno(h2Hom = 0.5, h2Het= [0,0], alpha = 0)
    est = h2Estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    return est 

@pytest.mark.C2S2_rand
def test_simNGCTA(C2S2_rand) :
    est = C2S2_rand
    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "GCTA", fixed_effects= ["Xc", "abcd_site"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

@pytest.mark.C2S2_rand
def test_simAdjHE(C2S2_rand):
    est = C2S2_rand
    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"], PC_effect = "fixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "AdjHE", fixed_effects= ["Xc"], random_groups = "abcd_site", PC_effect= "fixed")
    print(result)

    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"], PC_effect = "mixed")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "AdjHE", fixed_effects= ["Xc"], random_groups = "abcd_site", PC_effect= "mixed")
    print(result)

    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"], PC_effect = "random")
    print(result)
    result = est.estimate(mpheno = ["Y0"], npc = [0,1], Method = "AdjHE", fixed_effects= ["Xc"], random_groups = "abcd_site", PC_effect= "random")
    print(result)

    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05) 

#%%
@pytest.mark.C2S2_rand_thorough
def test_thoroughC2S2_rand(): 
    h2s = np.zeros((10,3,2))
    for i in range(10) : 
        sim = pheno_simulator(rng = i, nsubjects= 1000)
        sim.sim_sites(nsites =2, siteDistribution = "NOT_EQUAL")
        sim.sim_pops(nclusts= 2)
        sim.sim_genos()
        sim.sim_pheno(h2Hom = 0.5, h2Het= [0,0], alpha = 0)
        est = h2Estimation()
        est.GRM = sim.GRM
        est.df = sim.df
        h2s[i,0,0]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "AdjHE", random_groups="abcd_site", fixed_effects= ["Xc"], PC_effect = "fixed")["h2"][0]
        h2s[i,0,1]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"], PC_effect = "fixed")["h2"][0]
        h2s[i,1,0]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "AdjHE", random_groups="abcd_site", fixed_effects= ["Xc"], PC_effect = "mixed")["h2"][0]
        h2s[i,1,1]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"], PC_effect = "mixed")["h2"][0]
        h2s[i,2,0]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "AdjHE", random_groups="abcd_site", fixed_effects= ["Xc"], PC_effect = "random")["h2"][0]
        h2s[i,2,1]= est.estimate(mpheno = ["Y0"], npc = [1], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"], PC_effect = "random")["h2"][0]

