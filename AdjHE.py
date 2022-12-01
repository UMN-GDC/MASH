#! /usr/bin/env python3

"""
AdjHE estimator
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-06-06
"""

##############################################################
# The main file of the AdjHE estimator: loads, cleans, selects
# , loops, and store heritability estimates 
##############################################################

import os
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
#os.chdir("/panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit")
from functions.Data_input.parser import get_args, read_flags
from functions.traits_visualizer import covs_vs_cov_of_interest
from functions.Estimation.all_estimators import Basu_estimation

#######################################
#######################################
#%% For troubleshootingg
args= read_flags({"argfile" : "Example/basic_AdjHE.json"})
# args= read_flags({"argfile" : "Example/basic_GCTA.json"})
# args= read_flags({'argfile' : "simulations/Analysis_jsons/first_full.json"})
#######################################
#######################################


#%% Get command line arguments
# Get CL arguments and convert them to usable Python objects in a dictionary
args = read_flags(get_args())
print("These are the list of arguments that were input:")
print(args)

#%% Read in all data
ests = Basu_estimation(prefix = args["prefix"],
                          pheno_file = args["pheno"], 
                          cov_file= args["covar"], 
                          PC_file= args["PC"],
                          ids = args["ids"])
#%% Save images of covariate relations
# covs_vs_cov_of_interest(df, args["RV"], args["covars"], args["out"])

#%%
ests.looping(covars = args["covars"], npc = args["npc"], mpheno = args["mpheno"], loop_covars = args["loop_covars"])
ests.estimate(Method = "AdjHE", npc = [5])
# %%
print("Writing results")
ests.results.to_csv(args["out"], index=False, na_rep='NA')


