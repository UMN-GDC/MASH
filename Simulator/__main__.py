#! /usr/bin/env python3

"""
AdjHE estimator
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-06-06
"""

##############################################################
# Simulation runner for AdjHE package 
##############################################################

import os
from AdjHE.data_input.parser import get_args #, read_flags
from Simulator.simulation_helpers.run_sims import sim_experiment
import itertools
import json

#%% Get command line arguments
# Get CL arguments and convert them to usable Python objects in a dictionary
#args= read_flags({"argfile" : "Simulations/Sim_params/Adding_sites_clusts.json"})
args = get_args()
print("Simulating with the following parameters:")
print(args)


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
              races_differ=True,
              cov_effect=True,
              ortho_cov=True,
              random_BS=args["random_BS"])


# python program to check if a path exists
#if path doesn’t exist we create a new path
os.makedirs(os.path.dirname("Simulations/" + args["out"]+ ".csv"), exist_ok = True)

# write to out
df.dropna().to_csv(args["out"] + ".csv", 
          header=  True, index= False)


