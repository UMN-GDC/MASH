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
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
#os.chdir("/panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit")
import pandas as pd
import itertools
from functions.Data_input.load_data import load_everything
from functions.Data_input.parser import get_args, read_flags
from functions.traits_visualizer import covs_vs_cov_of_interest
from functions.Estimation.GCTA_wrapper import GCTA

#######################################
#######################################
#%% For troubleshootingg
# args= read_flags({"argfile" : "Example/basic_AdjHE.json"})
# args= read_flags({'argfile' : "simulations/Analysis_jsons/first_full.json"})
#######################################
#######################################


#%% Get command line arguments
# Get CL arguments and convert them to usable Python objects in a dictionary
args = read_flags(get_args())
print("These are the list of arguments that were input:")
print(args)


#%% Read in all data
df, GRM, phenotypes = load_everything(prefix = args["prefix"],
                          pheno_file = args["pheno"], 
                          cov_file= args["covar"], 
                          PC_file= args["PC"],
                          ids = args["ids"])
#%% Save images of covariate relations
# covs_vs_cov_of_interest(df, args["RV"], args["covars"], args["out"])

#%%
print("Calculating heritibility")

# Create list of covariate sets to regress over
if args["covars"] != None :
    # make them all lowercase
    args["covars"] = [covar.lower() for covar in args["covars"]]
    # Create the sets of covarates over which we can loop
    # This will return a list of lists of covariate names to regress on
    cov_combos = [args["covars"][0:idx+1] for idx, c in enumerate(args["covars"])]
    # If we don't want to loop, just grab the last item of the generated list assuming the user wants all of those variables included 
    if (args["loop_covs"] != True) : 
        cov_combos = [cov_combos[-1]]
else :
    cov_combos = [[]]

if args["mpheno"] == "all" :
    mpheno = phenotypes
else :
    # make them lowercase
    mpheno =  [ph.lower() for ph in args["mpheno"]]


#%%
# create empty list to store heritability estimates
results = pd.DataFrame()
covars = args["covars"]

if args["npc"] == None :
    args["npc"] = [0]

# Loop over each set of covariate combos
for covs, nnpc in itertools.product(cov_combos, args["npc"]):
    # Reset Precomputed values
    Precomputed = {}
    # For each set of covariates recalculate the projection matrix
    Precomputed = compute_Sjs(A, df, Precomputed, npc)
    Precomputed = compute_
    
    # loop over all combinations of pcs and phenotypes
    for mp, nnpc in itertools.product(mpheno, args["npc"]):
        ids = df[["fid", "iid"]]
        # seed empty result vector
        # result.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]
        # create the regression formula and columns for seelcting temporary
        form, cols  = create_formula(nnpc, covars, mp, RV)
        # save a temporary dataframe
        temp = df[cols].dropna()
        # Save residuals of selected phenotype after regressing out PCs and covars
        temp[mp] = smf.ols(formula = form, data = temp, missing = 'drop').fit().resid
        # Potentially could use this to control for random effects
        # smf.mixedlm(formula= form, data = temp, groups=temp["scan_site"])
        # keep portion of GRM without missingess for the phenotypes or covariates
        nonmissing = ids[ids.iid.isin(temp.iid)].index
        GRM_nonmissing = GRM[nonmissing,:][:,nonmissing]

        if not silent :
            print(Method + "Estimation...")
        # Select method of estimation
        if Method == "AdjHE": 
            result = load_n_AdjHE(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)

        elif Method == "MOM": 
            result = load_n_MOM(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)
        elif Method == "PredlMM" : 
            result = load_n_PredLMM(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)
        elif Method == "GCTA" :
            result = GCTA()

        
        r = load_n_estimate(
            df=df, covars=covs, nnpc=nnpc, mp=mp, GRM= GRM, std= False, Method = args["Method"], RV = args["RV"])
        results = pd.concat([results, r], ignore_index = True)
    
# %%
print("Writing results")
results.to_csv(args["out"], index=False, na_rep='NA')


