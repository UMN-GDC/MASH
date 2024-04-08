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
    parser = argparse.ArgumentParser(prog='Simulate data for genetic research.', description = """This program simulates genetic data for use in genetic research.
                                     while controlling for heritability and other factors. It is designed to be used with the AdjHE package.""")
    parser.add_argument('--argfile', help='json file with arguments to simulate data')
    parser.add_argument('--nSNPs', type=int, default=1000, help='Number of SNPs')
    parser.add_argument('--nsubjects', type=int, default=500, help='Number of subjects')
    parser.add_argument('--nclusts', type=int, default=1, help='Number of clusters')
    parser.add_argument('--nphenos', type=int, default=2, help='Number of phenotypes')
    parser.add_argument('--shared', type=float, default=0.5, help='Shared proportion')
    parser.add_argument('--prop_causal', nargs=2, type=float, default=[0.25, 0.25], help='Proportion of causal SNPs')
    parser.add_argument('--theta_alleles', nargs=2, type=float, default=[0.95, 0.25], help='Theta alleles')
    parser.add_argument('--h2Hom', type=float, default=0.8, help='Homogeneity')
    parser.add_argument('--h2Het', nargs=2, type=float, default=[0.1, 0.1], help='Heterogeneity')
    parser.add_argument('--savePlink', action='store_true', help='Save plink binary')
    parser.add_argument('--estimateHeritability', action='store_true', help='Estimate heritability') 
    args = vars(parser.parse_args())
    print(args)
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
    
    return(raw_args)






