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

def readable_file(path):
    """
    Throw argparse exception unless parameter is a valid readable filename 
    string. This is used instead of argparse.FileType("r") because the latter 
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
    except (AssertionError, OSError, TypeError) :
        raise argparse.ArgumentTypeError("Cannot read file at {}".format(path))
        
        
def readable_file_or_none(path):
    """
    Throw argparse exception unless path either is "none" or points to a 
    readable file.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename, or the string "none"
    """
    return path if path == "None" else readable_file(path)
        


