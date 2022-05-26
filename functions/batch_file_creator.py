#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AdjHE batch file maker
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-05-26
"""

##############################################################
# this creates batch files for parallel processing with 
# AdjHE estimator
##############################################################

# Read argfile
args={}
with open("Example/Batch_Arg_file.txt") as f:
    for line in f:
        (key, val) = line.split("=")
        # remove line break
        args[key] = val[:-1]
try:
    # convert a string of arugments sto a list
    args['mpheno'] = eval(args['mpheno'])
except:
    # Convert a single integer value to a list
    args['mpheno'] = list(args['mpheno'])





#%%

def single_batch(args, batch_out) :
    """
    Writes a batch file from the set of argsuments to the specified batch_output destination

    Parameters
    ----------
    args : dictionary
        a dictionary of AdjHE arguments.
    batch_out : str
        file path to save output batch file
        
    Returns
    -------
    None.

    """
    # open destination specified by out argument for writing 
    with open(batch_out, 'w') as f:
        # write file line by line separating arguments from name with =
        for key, value in args.items():
            f.write('%s=%s \n' % (key, value))

#%%
def split(a, n):
    """
    Split mphenos (a) into n approximately equal length lists

    Parameters
    ----------
    a : list
        list of integers for phenotypes.
    n : int
        number of batches to create.

    Returns
    -------
    list
        list of lists of phenotypes for batches.

    """
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))



# make all batches

def all_batches(args, nbatches) :
    """
    given an argfile and number of batches desired creates set of batch files

    Parameters
    ----------
    args : dict
        dictiponary of AdjHE arguments.
    nbatches : int
        number of desired batches to run.

    Returns
    -------
    None.

    """
    # split up phenotypes into nbatches roughly equally sized lists
    pheno_splits =  list(split(args["mpheno"], nbatches))
    
    # create temp args
    temp = args.copy()
    
    # loop over the desired number of batches
    for i in range(1, nbatches) :
        temp["mpheno"] = pheno_splits[i]
        single_batch(args = temp, batch_out= args["out"] +"_" + str(i) + "of" + str(nbatches))
        
        
    

#%%

all_batches(args, 9)


