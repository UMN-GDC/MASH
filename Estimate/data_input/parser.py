#! /usr/bin/env python3

"""
AdjHE command line argument parser
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-06-06
"""

##############################################################
# This is a helper tool for the AdjHE estimator that takes 
# command line arguments and converts them for use in AdjHE 
# heritability estimation
##############################################################

import argparse
import json 
from Estimate.data_input.types_n_valids import readable_file_or_none, readable_json

#%%

def get_args() :
    """
    Collects user input from command line and returns them for running with Estimate.

    Returns
    -------
    dictionary of user arguments from command line

    """
    # Create arg parser
    parser = argparse.ArgumentParser(
        prog='Running adjusted HE regression',description="This program gives estimation in formula fashion." 
        "Make sure you have enough memory to store GRM matrix in python.")
    
    
    # Required arguments: file paths for phenotypes, covariates, prinicpal components, GRM, and results output, 
    parser.add_argument(
        '--PC',
        type=readable_file_or_none,
        metavar= "EIGENVECTOR_FILE_PATH", 
        help='Read PLINK format covariate file contains the PCs'
        'PCs should be generated using the same set of individuals in GRM files.'
        'If --npc is not specified then all PCs in the file will be used.')
    
    parser.add_argument('--covar',
                        type=readable_file_or_none,
                        metavar= "COVARIATE_FILE_PATH", 
                        help='Read PLINK format covariate file contains covariates besides PCs to be adjusted')

    parser.add_argument('--prefix',
                        type=str,
                        metavar= "GRM_FILE_PREFIX", 
                        help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id')

    parser.add_argument('--pheno',
                        #type=argparse.FileType('r'),
                        type=readable_file_or_none,
                        metavar= "PHENOTYPE_FILE_PATH", 
                        help='Read PLINK format phenotype file [required]'
                        'If --mpheno is not specified then then 3rd column (the 1st phenotype) will be used.')

    parser.add_argument('--out',
                        type=str,
                        metavar= "OUTPUT_FILE_PATH", 
                        help='Specify the output file name.')

    # Optional
    parser.add_argument('--npc', 
                        nargs="+",
                        type=int,
                        metavar= "#_PCs", 
                        help='Specify the number of PCs to be adjusted for. Can be a list for model comparison')
    
    parser.add_argument('--mpheno',
                        nargs="+",
                        type=int,
                        metavar= "DESIRED_PHEN_INDEX", 
                        help='Specify which phenotype to use from phenotype file (Can be a list). The index starts after the FID and IID columns')
    
    parser.add_argument('--k',
                        type=int,
                        help='Specify the number of rows in restoring the GRM each time.'
                        'This could affect the computation time and memory especially when sample size is large. If not provide, it will process the whole GRM at one time.')
    
    parser.add_argument('--std',
                        action='store_true',
                        help='Run SAdj-HE (i.e., with standardization)')
    
    parser.add_argument('--covars',
                        nargs="+",
                        type=int,
                        metavar= "ORDERED_INDICES_OF_USEFUL_COVARS", 
                        help='Specify which covariates to control for from the covariate file. Should be a list of the column numbers not including the FID and IID columns')
    
    # Or accept a file with all arguments
    parser.add_argument("--argfile", 
                        # type=readable_json,
                        metavar= "ARGFILE_FILE_PATH", 
                        help="Filename to be passed containing all information for PC's, covariates, phenotypes, and grm")
    
    # Flag to loop over covariates or just do it once with all covariates
    parser.add_argument('--loop_covs',
                        action='store_true', 
                        help='Loop over the ordered list of covariates and retain all results.')
    
    # Flag to generate diagnostic plots w.r.t the assumptions required for the model
    parser.add_argument('--covar_relates',
                        action='store_true', 
                        help='Create and save diagnostic plots w.r.t the assumptions made in this model.')

    
    # Flag for using adjHE or PredLMM
    parser.add_argument('--RV', 
                        type=str,
                        help='Specify the random variable of interest')
    
    # Flag for using AdjHE simplification
    parser.add_argument('--fast', 
                        action='store_true',
                        help='Specify whether to use AdjHE estimation which is faster. Default is to use AdjHE method.')

    parser.add_argument('--Method',
                        type= str,
                        help = 'Specify which method of estimation you wish to use. (AdjHE, GCTA, PredLMM, SWD, COMBAT')
    
    # Set defaults
    parser.set_defaults(PC="None", fast = False, npc=None,
                        covar="None", mpheno=1, k=0,
                        prefix = "None", pheno = "None", out = "None",
                        PredLMM = False, RV = None, covar_relates = True,
                        loop_covs=False, argfile = None, covars =1,
                        std=  False)
    
    
    # return the arguments as a dictionary
    args = vars(parser.parse_args())
    # ForTroubleshooting  uncomment the next line
    # args['argfile'] = '/home/christian/Research/Stat_gen/tools/Basu_herit/Example/Arg_file.txt'
    return(args)




def read_flags(raw_args):
    """
    Takes the raw command line arguments and converts them to objects usable in Python.

    Parameters
    ----------
    raw_args : list
        raw arguments from the command line parser

    Returns
    -------
    dictionary of arguments usable in python for AdjHE

    """
    if raw_args['argfile'] != None :
        # Read data from standard json format
        with open(raw_args['argfile']) as f:
            raw_args = json.load(f)
    
    else : # Read each individual flag
        
        # Ensure types 
        raw_args["k"] = int(raw_args["k"])
        try:
            # convert a string of arugments sto a list
            raw_args['mpheno'] = eval(raw_args['mpheno'])
        except:
            # Convert a single integer value to a list
            raw_args['mpheno'] = list(raw_args['mpheno'])
          
        try:
            # convert a string of arugments sto a list
            raw_args['npc'] = eval(raw_args['npc'])
        except:
            # Convert a single integer value to a list
            raw_args['npc'] = list(raw_args['npc'])
            
        
            
        ## Do the same for specified covariates
        try:
            # Convert to    list of integers of agrument is a list
            raw_args["covars"] = eval(raw_args["covars"])
        except: 
            # convert single integer to integer list 
            raw_args["covars"] = raw_args["covars"]
    
    return(raw_args)



