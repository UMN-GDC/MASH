#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 19:06:54 2023

@author: christian
"""

from pathlib import Path
import pandas as pd
import numpy as np

#%%

path = 'Results/Het'  # or unix / linux / mac path

# Get the files from the path provided in the OP
files = Path(path).glob('**/*.csv')  # .rglob to get subdirectories

dfs = []
for file in files:
    if "Small_verifi" not in str(file):
        df = pd.read_csv(file)
        dfs.append(df)
        
pd.concat(dfs, ignore_index=True).to_csv(path + "/Het_IID_5000.csv",index= False, header= True, sep = ",")
