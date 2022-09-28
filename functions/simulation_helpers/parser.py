#! /usr/bin/env python3

"""
Simulation command line argument parser
Christian Coffman: coffm049@umn.edu
Created 2022-09-27
Last Updated 2022-09-27
"""

##############################################################
# This is a helper tool for the simualtions take 
# command line arguments 
##############################################################

import argparse
from functions.types_n_valids import readable_file_or_none

#%%

def get_sim_args() :
    """
    Collects user input from command line and returns them for simulating data.

    Returns
    -------
    dictionary of user arguments from command line

    """
    # Create arg parser
    parser = argparse.ArgumentParser(
        prog='Simulating data for heritability estimation',description="Simulating data for heritability estimation")
    
    
    # Required arguments: file paths for phenotypes, covariates, prinicpal components, GRM, and results output, 
    
    parser.add_argument('--covar',
                        type=readable_file_or_none,
                        metavar= "COVARIATE_FILE_PATH", 
                        help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')

    parser.add_argument('--prefix',
                        type=str,
                        metavar= "GRM_FILE_PREFIX", 
                        help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id')


    parser.add_argument('--out',
                        type=str,
                        metavar= "OUTPUT_FILE_PATH", 
                        help='Specify the output file name.')
    
    parser.add_argument("--sigmas", 
                        # type=readable_json,
                        metavar= "PATH_TO_SIGMA_LIST", 
                        help="Filename to be passed containing all sigma values to be simulated over.")    


    # Optional    
    parser.add_argument('--nreplicates',
                        type=int,
                        metavar= "NUM_REPS", 
                        help='Number of replicates at each parameter configuration')
            
    # Flag for using adjHE or PredLMM
    parser.add_argument('--RV', 
                        type=str,
                        help='Specify the random variable of interest')

    
    # Set defaults
    parser.set_defaults(covar=None, RV = None, argfile = None, covars =1)
    
    
    # return the arguments as a dictionary
    args = vars(parser.parse_args())
    return(args)







