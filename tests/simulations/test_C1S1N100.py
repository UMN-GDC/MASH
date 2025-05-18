#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 09:09:56 2023

@author: christian
"""
import itertools
import numpy as np
import statsmodels.api as sm
from Estimate.data_input.parser import get_args, read_flags
from Estimate.estimators.all_estimators import h2Estimation
from Simulate.simulation_helpers.run_sims import sim_experiment
import pytest

tolerance = 0.15

@pytest.fixture
def args() : 
    args = {
      "nsubjectss": [100],
      "sgs" : [0.25, 0.5, 0.75],
      "sss" : [0],
      "ses" : [0.25, 0.5, 0.75],
      "site_comps" : ["EQUAL"],
      "nsites" : [1],
      "theta_alleless" : [0.1],
      "nclustss" : [1],
      "dominances" : [3],
      "prop_causals" : [0.02],
      "site_deps" : [False],
      "nnpcs" :[0],
      "nSNPss" : [200],
      "phenss" : [1],
      "reps" : 5,
      "all_ests" : True,
      "site_het" : True,
      "Meta" : False,
      "site_2moment" : True,
      "random_BS" : True,
      "out" : "tests/simulations/results/S1C1N100_all"
    }

    return args 

@pytest.mark.C1S1N100
def test_C1S1N100(args) : 
    sigmas = []
    for sg, ss, se in itertools.product(args["sgs"], args["sss"], args["ses"]) :
        if sg + ss + se == 1 :
            if sg != 0 :
                sigmas += [[sg, ss, se]]
            elif (sg ==0) and (ss == 0) :
                sigmas += [[sg, ss, se]]
    
    
    # run experiments
    df = sim_experiment(nsubjectss = args["nsubjectss"],
                        sigmas = sigmas,
                        site_comps = args["site_comps"],
                        nsites = args["nsites"],
                        theta_alleless = args["theta_alleless"],
                        nclustss = args["nclustss"],
                        dominances= args["dominances"],
                        prop_causals= args["prop_causals"],
                        site_deps=args["site_deps"] ,
                        nnpcs = args["nnpcs"],
                        nSNPss= args["nSNPss"],
                        phenss= args["phenss"],
                        reps = args["reps"],
                        site_het = args["site_het"],
                        clusters_differ=True,
                        cov_effect=True,
                        ortho_cov=True,
                        random_BS=args["random_BS"])
    df = df.assign(sim_h2 = np.round(df["sg"]/(df["sg"] + df["se"]), decimals = 2))
    R2 = sm.OLS(df["sim_h2"], df["h2"]).fit().rsquared
    assert (1-tolerance) < R2
