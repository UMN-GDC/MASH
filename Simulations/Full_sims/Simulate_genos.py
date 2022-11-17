#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate phenotypes
Created on Thu Oct 13 09:25:44 2022

@author: christian
"""
from functions.simulation_helpers.Sim_generator import sim_experiment, plot_save
            
# %%
Sim = sim_experiment(nsubjectss = [1000], nsitess=[2], nclustss=[5], reps=3, nnpcs=[0])
