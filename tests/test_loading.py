#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 09:09:56 2023

@author: christian
"""
import os
import pytest
import numpy as np
import pandas as pd
from Estimate.data_input.parser import get_args, read_flags
from Estimate.data_input.load_data import ReadGRMBin
from Estimate.estimators.all_estimators import h2Estimation


       
        
#%%

# Call this with 
# python -m pytest -m loading 
@pytest.mark.loading
def test_loading_GRMBin() :
    n = 5
    l= np.tril_indices(n)

    df = np.random.uniform(-0.5, 0.5, size = len(l[0]))
    GRM1 = np.zeros((n,n))
    GRM1[l] = df
    GRM1 = GRM1 + GRM1.T - np.diag(np.diag(GRM1))
    
    GRM1[l].astype("f4").tofile("test.grm.bin")    

    pd.DataFrame({"a" : np.repeat(1, n), "b" : np.repeat(1,n)}).to_csv("test.grm.id", 
                                                                       sep = "\t", header= False, index = False)
    
    __, GRM2 = ReadGRMBin("test", sub_ids = None)
    
    np.testing.assert_allclose(GRM1, GRM2)
    os.remove("test.grm.bin")
    os.remove("test.grm.id")
    
@pytest.mark.loading
def test_loading_all() :
    args= read_flags({"argfile" : "AdjHE/examples/Generic.json"})
    ests = h2Estimation(prefix = args["prefix"],
                           pheno_file = args["pheno"], 
                           cov_file= args["covar"], 
                           PC_file= args["PC"],
                           ids = args["ids"])
    
    assert ['FID', 'IID', 'Head_size', 'Liver_size', 'O2_level', 'first100', 'first1000', 'first3000',
    'pc_1', 'pc_2', 'pc_3', 'pc_4', 'pc_5', 'pc_6', 'pc_7', 'pc_8', 'pc_9', 'pc_10', 'Tidyness', 
    'Happy_camper', 'Like_of_levis', 'abcd_site'] == ests.df.columns.tolist()


@pytest.mark.loading
def test_argfile_iid_fid_cols():
    """Test that iid_col and fid_col are correctly read from argfile."""
    args = read_flags({"argfile": "tests/test_data/config_AdjHE.json"})
    assert args.get("iid_col") == "IID"
    assert args.get("fid_col") == "FID"


@pytest.mark.loading
def test_argfile_default_iid_fid():
    """Test default iid_col and fid_col when not specified in argfile."""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--iid-col', type=str, default='IID')
    parser.add_argument('--fid-col', type=str, default='FID')
    args = parser.parse_args([])
    assert args.iid_col == 'IID'
    assert args.fid_col == 'FID'


@pytest.mark.loading
def test_argfile_custom_iid_fid_cols():
    """Test that custom iid_col and fid_col can be specified in argfile."""
    custom_config = {
        "PC": "tests/test_data/EUR2.2.eigenvec",
        "covar": ["tests/test_data/EUR_covariates.cov"],
        "prefix": "tests/test_data/EUR2",
        "pheno": "tests/test_data/EUR_simulation2.pheno",
        "out": "/users/4/coffm049/software/MASH/tests/results/EUR_custom_cols",
        "npc": [2],
        "mpheno": [1],
        "k": 0,
        "std": False,
        "argfile": None,
        "loop_covars": False,
        "covar_relates": False,
        "random_groups": None,
        "Naive": False,
        "Method": "AdjHE",
        "ids": None,
        "qcovar": None,
        "covar_discrete": None,
        "iid_col": "SubjectID",
        "fid_col": "FamilyID"
    }
    import json
    with open("/users/4/coffm049/software/MASH/tests/test_data/config_custom_ids.json", "w") as f:
        json.dump(custom_config, f)
    
    args = read_flags({"argfile": "/users/4/coffm049/software/MASH/tests/test_data/config_custom_ids.json"})
    assert args.get("iid_col") == "SubjectID"
    assert args.get("fid_col") == "FamilyID"

    

        
    
    
        





