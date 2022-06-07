#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 08:25:47 2022

@author: christian
"""

#TODO: 
#   - check for GRM files readable
#   - check phenos, covars, pcs having 3+ columns


import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
covar = "/home/christian/Research/Stat_gen/tools/Basu_herit/Example/covar.txt"

#%%

def three_col_file(path):
    """
    Throw argparse exception unless parameter is a valid readable filename 
    string with more than three columns. This is used instead of argparse.FileType("r") because the latter 
    leaves an open file handle, which has caused problems.
    :param 
    :return: 
    Parameters
    ----------
    path: str
        Parameter to check if it represents a valid filename.

    Returns
    -------
    str String representing a valid filename

    """
    try :
        assert os.access(path, os.R_OK)
        return os.path.abspath(path)
        df = pd.read_table(path, sep = " ", header=0,  nrow=3) 
        # if fewer than 3 columns raise an error
        if df.shape[1] < 3 :
            raise argparse.ArgumentTypeError(
                "File at {} has fewer than three columns. Are you sure this file contians FID, IID, and then hte rest of the data?".format(path))
        else :
            return None
    except (AssertionError, OSError, TypeError) :
        raise argparse.ArgumentTypeError("Cannot read file at {}".format(path))
#%%        


def readable_file_or_none(path):
    """
    Throw argparse exception unless path either is "none" or points to a 
    readable file.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename, or the string "none"
    """
    return path if path == "None" else three_col_file(path)




