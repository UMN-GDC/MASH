from AdjHE.estimation.GCTA_wrapper import GCTA, gcta
from Simulator.simulation_helpers.Sim_generator import pheno_simulator
import pytest
#%%

# Call this with 
# python -m pytest -m sim_working 
@pytest.mark.sim_working
def test_loading_all() :
    sim = pheno_simulator(nsubjects= 500)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 1)
    sim.sim_genos()
    sim.sim_gen_effects()
    sim.sim_covars()
    sim.sim_pheno(var_comps = [0.5,0.0, 0.5])
    result = GCTA(
        df = sim.df,
        covars = ["Xc"],
        nnpc = 1 ,
        mp = "Y0",
        GRM = sim.GRM,
        gcta  = gcta,
        silent = False)
    assert result["h2"][0] == 0.5


