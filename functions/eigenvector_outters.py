#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:57:20 2022

@author: christian
"""

import numpy as np
import pandas as pd
import regex as re
# #%%
# subjects = 10
# nnpc = 5
# other_cols = 4

# X =  pd.DataFrame(np.random.rand(subjects, nnpc + other_cols))
# X.columns = ["pc_" + str(i) for i in range(1, nnpc + 1)] + ["a", "b", "c", "d"]

# #%%
# # seed the empty lambda matrix 
# Lambda_mat = pd.DataFrame(np.zeros((subjects ** 2, 0)))

# #%%

# r = re.compile("pc_*")
# newlist = list(filter(r.match, X.columns)) # Read Note below
# #%%
# df = pd.DataFrame(np.zeros(( int((subjects**2+ subjects)/2), 0)))

# # outer product for each individual column
# for c_ind in range(1, nnpc + 1) :
#     # construct the pc name
#     pc_col = "pc_" + str(c_ind)
#     t= np.outer(X[pc_col], X[pc_col])[np.triu_indices(subjects)]
#     df[str(c_ind)] = t

# #%%


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

