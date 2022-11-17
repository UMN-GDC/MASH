#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
from functions.simulation_helpers.Sim_generator import sim_experiment
            
# %%
Sim = sim_experiment(nsubjectss = [500, 1000, 5000], nsitess=[30], nclustss=[2,5], reps=50, nnpcs=[2])


Sim.to_csv("Simulations/Full_sims/Full_sim.csv", header= True, index= False)
