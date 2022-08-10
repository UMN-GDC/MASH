#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 08:25:47 2022

@author: christian
"""

#TODO: 
#   - check for GRM files readable


import pandas as pd
import argparse
import os

#parser = argparse.ArgumentParser()
#covar = "/home/christian/Research/Stat_gen/tools/Basu_herit/Example/covar.txt"

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
    except (AssertionError, OSError, TypeError) :
        raise argparse.ArgumentTypeError("Cannot read file at {}".format(path))
        
    # if fewer than 3 columns raise an error
    try:
        df = pd.read_table(path, sep = " ", header=0,  nrow=3) 
        assert df.shape[1] > 2
        return os.path.abspath(path)
    except(AssertionError, TypeError) :
        raise argparse.ArgumentTypeError(
            "File at {} has fewer than three columns. Are you sure this file contians FID, IID, and then hte rest of the data?".format(path))

#%%        


def readable_file_or_none(path):
    """
    Throw argparse exception unless path either is "none" or points to a 
    readable file.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename, or the string "none"
    """
    return path if path == "None" else three_col_file(path)



def readable_json(path):
    """
    Throw argparse exception unless path either is "none" or points to a 
    readable file.
    :param path: Parameter to check if it represents a valid filename
    :return: String representing a valid filename, or the string "none"
    """

    try :
        assert os.access(path, os.R_OK)
    except (AssertionError, OSError, TypeError) :
        raise argparse.ArgumentTypeError("Cannot read file at {}".format(path))

    try :
        base, ext =  os.path.splitext(path)
        assert ext == ".json"
    except (AssertionError, OSError, TypeError) :
        raise argparse.ArgumentTypeError("File at {} is not a .json file".format(path))

 
