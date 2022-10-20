#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:57:20 2022

@author: christian
"""

import numpy as np
import pandas as pd
import regex as re
#%%

                          
def multiple_outer(df, nnpc) :
    """
    Function for calculating multiple outer products for the case when more than one eigenvector is used.

    Parameters
    ----------
    df : pandas dataframe 
        n x p contianing FID, IID, and PC's

    Returns
    -------
    pandas dataframe
        n**2 x p
        each column represents the flattened upper half of the outer product of each PC.

    """
    subjects = df.shape[0]
    
    # get list of columns that are PCs
    # r = re.compile("pc_*")
    # pc_cols = list(filter(r.match, df.columns))
    
    # seed empty dataframe for outer product
    outers = pd.DataFrame(np.zeros(( int((subjects**2+ subjects)/2), 0)))

    # Get outter product of each PC column
    for c_ind in range(1, nnpc + 1) :
        # construct the pc name
        pc_col = "pc_" + str(c_ind)
        t= np.outer(df[pc_col], df[pc_col])[np.triu_indices(subjects)]
        outers["pc_" + str(c_ind)] = t
    # return outter product matrix
    return(outers)


