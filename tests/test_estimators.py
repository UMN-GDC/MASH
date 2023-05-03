#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 09:09:56 2023

@author: christian
"""
from AdjHE.data_input.parser import get_args, read_flags
from AdjHE.estimation.all_estimators import Basu_estimation


import pytest


# all pytests will read the same json file
@pytest.fixture
def args() : 
    return read_flags({"argfile" : "AdjHE/examples/Generic.json"})

# all pytests have the common fixture of reading arguments
@pytest.fixture
def estimator(args) :
    ests = Basu_estimation(prefix = args["prefix"],
                           pheno_file = args["pheno"],
                           cov_file = args["covar"],
                           PC_file = args["PC"],
                           ids = args["ids"])
    return ests


#%%

# Call this with 
# python -m pytest -m loading 
@pytest.mark.AdjHE
def test_AdjHE_basic(args, estimator) : 
    estimator.estimate(Method = "AdjHE", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"], 
                  random_groups = args["random_groups"], Naive= args["Naive"])
    h = estimator.results["h2"][0]
    
    assert (h<0.9) and (h > 0.5)

@pytest.mark.AdjHE
def test_AdjHE_w_site_effects(args, estimator) :
    estimator.estimate(Method = args["Method"], npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= args["Naive"])
    h= estimator.results["h2"][0]
    assert (h<0.9) and (h > 0.5)

@pytest.mark.AdjHE
def test_AdjHE_naive(args, estimator) :
    estimator.estimate(Method = args["Method"], npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= True)
    h = estimator.results["h2"][0]
    assert (h<0.9) and (h>0.5)




