import pytest
import numpy as np
import pandas as pd
from Simulator.simulation_helpers.sim_effects import sim_effects

tolerance = 0.01

@pytest.fixture
def rng() :
    rng = np.random.default_rng(123)
    return rng


@pytest.fixture
def Genotypes(rng): 
    X = rng.binomial(2,0.5, size = (100, 50))
    return X

@pytest.mark.sim_effects
@pytest.mark.parametrize("freqDep,space", [(False, False),
                                                  (False, True),
                                                  (True, False)])
def test_sim_effects(rng, Genotypes, freqDep, space) :
    target = 0.5
    contrib = sim_effects(rng, Genotypes, target, 
                          freqDep, space) 

    assert abs(np.var(contrib) - target) < tolerance

