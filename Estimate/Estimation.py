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
# args = read_flags({"argfile": "/panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/03_Herit_ests/Asegs/full/New/Covbat.json"})
import os
import logging
import pandas as pd
#os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
# os.chdir("/panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit")
from Estimate.data_input.parser import get_args, read_flags
# from Estimate.traits_visualizer import covs_vs_cov_of_interest
from Estimate.estimators.all_estimators import h2Estimation

def main():
    #%% Get command line arguments
    # Get CL arguments and convert them to usable Python objects in a dictionary
    #args = read_flags({"argfile" : "../GCTA.json"})
    args = read_flags(get_args())
    
    logging.basicConfig(filename=f'{args["out"]}.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s',
                        level=logging.DEBUG)
    logging.info("Started")
    logging.info("These are the list of arguments that were input:")
    logging.info(args)
    logging.info("Loading data")
    #%% Read in all data
    
    ests = h2Estimation(args = args)
    #%% Save images of covariate relations
    # covs_vs_cov_of_interest(df, args["RV"], args["covars"], args["out"])
    
    #%%
    logging.info("Estimating")
    ests.df = pd.merge(ests.ids, ests.df, on = ["FID", "IID"], how = "left")
    ests.estimate()
    # %%
    logging.info(f"Writing results to {args['out']}")
    ests.results.to_csv(args["out"], index=False, na_rep='NA')
    logging.info("Finished")
    
