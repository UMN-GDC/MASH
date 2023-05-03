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
import logging
from AdjHE.data_input.parser import get_args, read_flags
# from AdjHE.traits_visualizer import covs_vs_cov_of_interest
from AdjHE.estimation.all_estimators import Basu_estimation

#######################################
#######################################
#%% For troubleshootingg
# args= read_flags({"argfile" : "Example/basic_AdjHE.json"})
# args= read_flags({"argfile" : "Example/basic_GCTA.json"})
# args= read_flags({'argfile' : "simulations/Analysis_jsons/first_full.json"})
#######################################
#######################################


if __name__ == "__main__" :
    #%% Get command line arguments
    # Get CL arguments and convert them to usable Python objects in a dictionary
    args = read_flags(get_args())
    
    logging.basicConfig(filename=f'{args["out"]}.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s',
                        level=logging.DEBUG)
    logging.info("Started")
    logging.info("These are the list of arguments that were input:")
    logging.info(args)
    logging.info("Loading data")
    #%% Read in all data
    
    ests = Basu_estimation(prefix = args["prefix"],
                              pheno_file = args["pheno"], 
                              cov_file= args["covar"], 
                              PC_file= args["PC"],
                              ids = args["ids"])
    #%% Save images of covariate relations
    # covs_vs_cov_of_interest(df, args["RV"], args["covars"], args["out"])
    
    #%%
    logging.info("Estimating")
    ests.estimate(Method = args["Method"], npc = args["npc"], fixed_effects = args["fixed_effects"], mpheno = args["mpheno"], loop_covars = args["loop_covars"], random_groups = args["random_groups"], Naive= args["Naive"])
    # %%
    logging.info(f"Writing results to {args['out']}")
    ests.results.to_csv(args["out"], index=False, na_rep='NA')
    logging.info("Finished")
    
