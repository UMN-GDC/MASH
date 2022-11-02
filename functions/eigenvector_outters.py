#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 11:57:20 2022

@author: christian
"""

import numpy as np
                          
def columnwise_outter(mat) :
    """
    Gvien a (n x p) matrix of eigenvector loadings for individuals with p eigen loading vectors as the columns and n observations
    calculate and return all elementwise outer products stacked.

    Parameters
    ----------
    mat : numpy matrix
        (n x p) matrix with p columns for each eigenvetor loading and n observations.

    Returns
    -------
    St : np.array
        (p x n x n) stacked array of eigenvector loading outer products.

    """
    St = np.zeros((mat.shape[1], mat.shape[0], mat.shape[0]))

    for i in range(mat.shape[1]) :
        St[i,:,:] = np.outer(mat[:,i], mat[:,i])
    return St


