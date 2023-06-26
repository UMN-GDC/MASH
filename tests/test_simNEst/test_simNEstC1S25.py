from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from AdjHE.estimation.all_estimators import Basu_estimation
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import pytest

#%% Testing for Multiple sites and single cluster
@pytest.fixture
def C1S25() :
    rng = np.random.default_rng(12)
    sim = pheno_simulator(rng = rng, nsubjects= 1000)
    sim.sim_sites(nsites =25)
    sim.sim_pops(nclusts= 1)
    sim.sim_genos()
    sim.sim_gen_effects()
    sim.sim_covars()
    sim.sim_pheno(h2 = 0.5)
    est = Basu_estimation()
    est.GRM = sim.GRM
    est.df = sim.df
    return est

@pytest.mark.C1S25
def test_simNGCTA(C1S25) :
    est = C1S25
    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "GCTA", fixed_effects= ["Xc", "abcd_site"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05)

@pytest.mark.C1S25
def test_simNAdjHE(C1S25) :
    est = C1S25
    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "AdjHE", fixed_effects= ["Xc", "abcd_site"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05)

@pytest.mark.C1S25
def test_simNAdjHERE(C1S25):
    est= C1S25
    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "AdjHE", random_groups="abcd_site", fixed_effects= ["Xc"])
    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05)


#%% Testin for alternative estimators
#@pytest.mark.C1S25
#def test_simNSWD(C1S25):
#    est= C1S25
#    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "SWD", random_groups="abcd_site", fixed_effects= ["Xc"])
#    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05)

#@pytest.mark.C1S25
#def test_simNCombat(C1S25):
#    est= C1S25
#    result = est.estimate(mpheno = ["Y0"], npc = [0], Method = "COMBAT", random_groups="abcd_site", fixed_effects= ["Xc"])
#    assert result["h2"][0] == pytest.approx(0.5, abs = 0.05)

