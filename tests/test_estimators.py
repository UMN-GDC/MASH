#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 09:09:56 2023

@author: christian
"""
from AdjHE.data_input.parser import get_args, read_flags
from AdjHE.estimation.all_estimators import Basu_estimation


import pytest


@pytest.fixture
def args() : 
    return read_flags({"argfile" : "AdjHE/examples/Generic.json"})

@pytest.fixture
def estimator(args) :
    ests = Basu_estimation(prefix = args["prefix"],
                           pheno_file = args["pheno"],
                           cov_file = args["covar"],
                           PC_file = args["PC"],
                           ids = args["ids"])
    return ests


@pytest.mark.AdjHE
def test_AdjHE_basic(args, estimator) : 
    estimator.estimate(Method = "AdjHE", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"], 
                  random_groups = None, Naive= False)
    h = estimator.results["h2"][0]
    
    assert (h<0.9) and (h > 0.5)

@pytest.mark.AdjHE
def test_AdjHE_RE(args, estimator) :
    estimator.estimate(Method = "AdjHE", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= False)
    h = estimator.results["h2"][0]
    
    assert (h<0.9) and (h > 0.5)

@pytest.mark.AdjHE
def test_AdjHE_meta(args, estimator) :
    estimator.estimate(Method = "AdjHE", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= True)
    h = estimator.results["h2"][0]
    
    assert (h<0.9) and (h > 0.5)

# Test GCTA

@pytest.mark.GCTA
def test_GCTA_basic(args, estimator) :
    estimator.estimate(Method = "GCTA", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= False)
    h = estimator.results["h2"][0]
    assert (h<0.9) and (h>0.5)

@pytest.mark.GCTA
def test_GCTA_meta(args, estimator) :
    estimator.estimate(Method = "GCTA", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= True)
    h = estimator.results["h2"][0]
    assert (h<0.9) and (h>0.5)


# Test SWD
@pytest.mark.SWD
def test_SWD(args, estimator) :
    estimator.estimate(Method = "SWD", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= False)
    h = estimator.results["h2"][0]
    assert (h<0.9) and (h>0.5)

# Test Combat
@pytest.mark.Combat
def test_Combat(args, estimator) :
    estimator.estimate(Method = "Combat", npc = [3], fixed_effects = args["fixed_effects"],
                  mpheno = args["mpheno"], loop_covars = args["loop_covars"],
                  random_groups = "abcd_site", Naive= False)
    h = estimator.results["h2"][0]
    assert (h<0.9) and (h>0.5)




