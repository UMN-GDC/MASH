#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 08:50:21 2022

@author: christian
"""

import os
import numpy as np
os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit")
#os.chdir("/panfs/roc/groups/3/rando149/coffm049/tools/Basu_herit")
import pandas as pd
from functions.load_data import load_everything
from functions.parser import read_flags
from functions.eigenvector_outters import columnwise_outter
import statsmodels.formula.api as smf


#%% Define AdjHE estimator
