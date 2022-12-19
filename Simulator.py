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
#os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
from functions.Data_input.sim_parser import get_args
from run_sims import sim_experiment
import itertools


#%% Get command line arguments
# Get CL arguments and convert them to usable Python objects in a dictionary
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
	      all_ests = args["all_ests"])


# python program to check if a path exists
#if path doesn’t exist we create a new path
os.makedirs(os.path.dirname("Simulations/" + args["out"]+ ".csv"), exist_ok = True)

# write to out
df.dropna().to_csv("Simulations/" + args["out"] + ".csv", 
          header=  True, index= False)


