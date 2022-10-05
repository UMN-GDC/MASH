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
#os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
import pandas as pd
import itertools
from functions.AdjHE_estimator import load_n_estimate
from functions.load_data import load_everything
from functions.parser import get_args, read_flags
from functions.traits_visualizer import covs_vs_cov_of_interest

#######################################
#######################################
#%% For troubleshootingg
#c_args= {}
# c_args['argfile'] = "Example/Argfile.json"
#c_args['argfile'] = "simulations/sims_Argfile.json"
#args = read_flags(c_args)
#######################################
#######################################


#%% Get command line arguments
# Get CL arguments and convert them to usable Python objects in a dictionary
args = read_flags(get_args())
print("These are the list of arguments that were input:")
print(args)


#%% Read in all data
(df, covariates, phenotypes, GRM_array_nona, ids) = load_everything(prefix = args["prefix"],
                                                                    pheno_file = args["pheno"], 
                                                                    cov_file= args["covar"], 
                                                                    PC_file= args["PC"],
                                                                    k= args["k"])
#%% Save images of covariate relations
# covs_vs_cov_of_interest(df, args["RV"], args["covars"], args["out"])

#%%
print("Calculating heritibility")

# Grab covariate names if covariates specified
if (args["covars"] != None) and (len(covariates) > 0) :
    covars = [covar-1 for covar in args["covars"]]
    # Create the sets of covarates over which we can loop
    # This will return a list of lists of indices for the sets of covaraites to use
    cov_combos = [covars[0:idx+1] for idx, c in enumerate(covars)]
    # This will create a list of lists of the actually covariate names to control for
    cov_combos = [list(covariates[cov_combo]) for cov_combo in cov_combos]
    # If we don't want to loop, just grab the last item of the generated list assuming the user wants all of those variables included 
    if (args["loop_covs"] != True) : 
        cov_combos = [cov_combos[-1]]

else :
    covars = None

if isinstance(args["mpheno"], str) :
    mpheno = phenotypes
else :

    #%% get list of phenotype names to regress
    if (args["pheno"] != None) :
        # if separate phenotype file is specified, grab from that
        mpheno = [phenotypes[i-1] for i in args["mpheno"]] 
    else :
        # if separate pheno file is not specified grab from the covariates file
        mpheno = [covariates[i-1] for i in args["mpheno"]]


#%%
# create empty list to store heritability estimates
results = pd.DataFrame()

# loop over all combinations of pcs and phenotypes
if (covars == None) and (args["npc"] != None):
    for mp, nnpc in itertools.product(mpheno, args["npc"]):
        r = load_n_estimate(
            df=df, covars=[], nnpc=nnpc, mp=mp, ids=ids, GRM_array_nona=GRM_array_nona, std= False, fast = args["fast"])
        results = pd.concat([results, r])
elif (covars == None) and (args["npc"] == None):
    for mp in mpheno:
        r = load_n_estimate(
            df=df, covars=[], nnpc=0, mp=mp, ids=ids, GRM_array_nona=GRM_array_nona, std= False, fast = args["fast"])
        results = pd.concat([results, r])
elif (covars != None) and (args["npc"] == None):
    for mp, covs in itertools.product(mpheno, cov_combos):
        r = load_n_estimate(
            df=df, covars= covs, nnpc=0, mp=mp, ids=ids, GRM_array_nona=GRM_array_nona, std= False, fast = args["fast"])
        results = pd.concat([results, r])
else:
    
    for mp, nnpc, covs in itertools.product(mpheno, args["npc"], cov_combos):
        r = load_n_estimate(
            df=df, covars=covs, nnpc=nnpc, mp=mp, ids=ids, GRM_array_nona=GRM_array_nona, std= False, fast = args["fast"])
        results = pd.concat([results, r])
    

# %%
print("Writing results")
results.to_csv(args["out"] + "/Results/AdjHE.csv", index=False, na_rep='NA')
