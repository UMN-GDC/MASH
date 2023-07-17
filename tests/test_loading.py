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
from Estimate.estimators.all_estimators import Basu_estimation


       
        
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
    ests = Basu_estimation(prefix = args["prefix"],
                           pheno_file = args["pheno"], 
                           cov_file= args["covar"], 
                           PC_file= args["PC"],
                           ids = args["ids"])
    
    assert ['FID', 'IID', 'Head_size', 'Liver_size', 'O2_level', 'first100', 'first1000', 'first3000',
    'pc_1', 'pc_2', 'pc_3', 'pc_4', 'pc_5', 'pc_6', 'pc_7', 'pc_8', 'pc_9', 'pc_10', 'Tidyness', 
    'Happy_camper', 'Like_of_levis', 'abcd_site'] == ests.df.columns.tolist()

    

        
    
    
        





