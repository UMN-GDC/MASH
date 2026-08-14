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

def mpheno_type(value):
    """
    Argparse type for --mpheno: accept integer phenotype indices or the
    string "ALL" (case-insensitive) to run across all phenotypes.
    """
    if str(value).lower() == "all":
        return "ALL"
    return int(value)

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
                        help='Read covariate file(s). Supports comma-separated list for multiple files. '
                        'Format detected by extension: .csv (comma), .tsv/.tab (tab), whitespace (PLINK). '
                        'Column headers are preserved and can be used in --covars, --qcovar, --covar_discrete.')

    parser.add_argument('--prefix',
                        type=str,
                        metavar= "GRM_FILE_PREFIX", 
                        help='prefix for GCTA format GRM files, including PREFIX.grm.bin, PREFIX.grm.N.bin, and PREFIX.grm.id')

    parser.add_argument('--pheno',
                        #type=argparse.FileType('r'),
                        type=readable_file_or_none,
                        metavar= "PHENOTYPE_FILE_PATH", 
                        help='Read phenotype file [required]. Format detected by extension: '
                        '.csv (comma), .tsv/.tab (tab), .parquet, or whitespace (PLINK). '
                        'First two columns should be FID and IID, followed by phenotype columns. '
                        'Column headers recommended but not required.')

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
                        type=mpheno_type,
                        metavar= "DESIRED_PHEN_INDEX", 
                        help='Specify which phenotype to use from phenotype file (Can be a list). The index starts after the FID and IID columns. '
                        'Use "ALL" to run across all phenotypes.')
    
    parser.add_argument('--k',
                        type=int,
                        help='Specify the number of rows in restoring the GRM each time.'
                        'This could affect the computation time and memory especially when sample size is large. If not provide, it will process the whole GRM at one time.')
    
    parser.add_argument('--std',
                        action='store_true',
                        help='Run SAdj-HE (i.e., with standardization)')
    
    parser.add_argument('--covars',
                        nargs="+",
                        type=str,
                        metavar= "COVARIATE_NAMES", 
                        help='Space-separated list of covariate column names to adjust for')
    
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
    
    parser.add_argument('--qcovar',
                        nargs="+",
                        type=str,
                        metavar= "QCOVAR_LIST", 
                        help='Space-separated list of quantitative covariates for GCTA')
    
    parser.add_argument('--covar_discrete',
                        nargs="+",
                        type=str,
                        metavar= "DISCRETE_COVAR_LIST", 
                        help='Space-separated list of discrete covariates for GCTA')
    
    parser.add_argument('--pheno-filter',
                        type=str,
                        metavar= "PHENO_FILTER",
                        help='Filter phenotype data by condition (e.g., "age>30", "sex==1"). '
                        'Supports: ==, !=, >, <, >=, <= and column names from phenotype file.')
    
    parser.add_argument('--covar-filter',
                        type=str,
                        metavar= "COVAR_FILTER",
                        help='Filter by covariate values (e.g., "sex==1", "age>=18"). '
                        'Supports: ==, !=, >, <, >=, <= and column names from covariate file.')

    parser.add_argument('--iid-col',
                        type=str,
                        default='IID',
                        metavar= "IID_COLUMN",
                        help='Name of the IID column in phenotype/covariate files. Default: IID')

    parser.add_argument('--fid-col',
                        type=str,
                        default='FID',
                        metavar= "FID_COLUMN",
                        help='Name of the FID column in phenotype/covariate files. Default: FID')

    parser.add_argument('--na-values',
                        nargs="+",
                        metavar= "NA_VALUE",
                        help='Additional values to recognize as NA/missing in delimited input files '
                        '(e.g., -777 -888). Accepts a list of values.')

    # Set defaults
    parser.set_defaults(PC="None", fast = False, npc=None,
                        covar="None", mpheno=1, k=0,
                        prefix = "None", pheno = "None", out = "None",
                        PredLMM = False, RV = None, covar_relates = True,
                        loop_covs=False, argfile = None, covars =1,
                        std= False, qcovar=None, covar_discrete=None,
                        pheno_filter=None, covar_filter=None,
                        iid_col='IID', fid_col='FID',
                        na_values=None)
    
    
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
            raw_args['npc'] = eval(raw_args['npc'])
        except:
            # Convert a single integer value to a list
            raw_args['npc'] = list(raw_args['npc'])
        
        # Handle qcovar and covar_discrete - convert to list or None
        for key in ['qcovar', 'covar_discrete']:
            if key in raw_args and raw_args[key] is not None:
                if isinstance(raw_args[key], str):
                    raw_args[key] = raw_args[key].split(',')
                elif isinstance(raw_args[key], list):
                    pass  # Already a list
                else:
                    raw_args[key] = None
            
        ## Do the same for specified covariates
        # try:
        #     # Convert to    list of integers of agrument is a list
        #     raw_args["covars"] = eval(raw_args["covars"])
        # except: 
        #     # convert single integer to integer list 
        #     raw_args["covars"] = raw_args["covars"]
    
    # Normalize mpheno ("ALL" or list of ints) across both config and CLI paths
    raw_args['mpheno'] = normalize_mpheno(raw_args.get('mpheno'))
    # Normalize na_values to a list of scalars (or None)
    raw_args['na_values'] = normalize_na_values(raw_args.get('na_values'))
    
    return(raw_args)


def normalize_mpheno(mpheno):
    """
    Ensure mpheno is either the string "ALL" (case-insensitive), a list of
    integer phenotype indices, or a list of phenotype column names. Applies to
    both config (JSON argfile) and CLI.
    """
    if mpheno is None:
        return None
    if isinstance(mpheno, str):
        if mpheno.lower() == "all":
            return "ALL"
        try:
            return [int(mpheno)]
        except ValueError:
            return [mpheno]  # a single column name
    if isinstance(mpheno, list):
        if len(mpheno) == 1 and isinstance(mpheno[0], str) and mpheno[0].lower() == "all":
            return "ALL"
        out = []
        for m in mpheno:
            if isinstance(m, str):
                try:
                    out.append(int(m))
                except ValueError:
                    out.append(m)  # keep column names as-is
            else:
                out.append(m)
        return out
    return [int(mpheno)]  # single scalar


def normalize_na_values(na_values):
    """
    Ensure na_values is a list of scalars (or None) so it can be passed
    directly to pandas read functions. Numeric-looking strings (e.g. from
    the command line) are converted to int/float so they match numeric
    columns during parsing.
    """
    if na_values is None:
        return None
    if not isinstance(na_values, list):
        na_values = [na_values]
    out = []
    for v in na_values:
        if isinstance(v, str):
            try:
                out.append(int(v))
            except ValueError:
                try:
                    out.append(float(v))
                except ValueError:
                    out.append(v)
        else:
            out.append(v)
    return out






