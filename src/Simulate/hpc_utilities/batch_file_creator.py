#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AdjHE batch file maker
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-05-26
"""

##### TODO 
# Still need to decide where it will put these files after making them, and then if it will delete them after use in a slurm script

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
    args['continuousPhenos'] = eval(args['continuousPhenos'])
except:
    # Convert a single integer value to a list
    args['continuousPhenos'] = list(args['continuousPhenos'])

try:
    args['binPhenos'] = eval(args.get('binPhenos', 'None'))
except:
    args['binPhenos'] = None

try:
    args['prevalence'] = eval(args.get('prevalence', 'None'))
except:
    args['prevalence'] = None

# Ensure preprocess is set
if 'preprocess' not in args:
    args['preprocess'] = 'None'




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
    Split phenotypes (a) into n approximately equal length lists

    Parameters
    ----------
    a : list
        list of phenotype indices or names.
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
        dictionary of AdjHE arguments.
    nbatches : int
        number of desired batches to run.

    Returns
    -------
    None.

    """
    # split up continuous phenotypes into nbatches roughly equally sized lists
    continuous_splits = list(split(args["continuousPhenos"], nbatches))
    
    # split up binary phenotypes into nbatches roughly equally sized lists
    bin_splits = list(split(args.get("binPhenos", []), nbatches)) if args.get("binPhenos") else []
    
    # create temp args
    temp = args.copy()
    
    # loop over the desired number of batches
    for i in range(0, nbatches) :
        if i < len(continuous_splits):
            temp["continuousPhenos"] = continuous_splits[i]
        else:
            temp["continuousPhenos"] = []
        if i < len(bin_splits):
            temp["binPhenos"] = bin_splits[i]
        else:
            temp["binPhenos"] = []
        single_batch(args = temp, batch_out= args["out"] +"_" + str(i + 1) + "of" + str(nbatches))
