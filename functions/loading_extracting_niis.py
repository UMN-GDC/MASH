#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 13:20:31 2022

@author: christian
"""
import os
import numpy as np
import nibabel as nib
import pandas as pd
import re
#%%

def load_extract_niis(files, save = False, save_path = None) :
    files = pd.read_csv(files, header=None)

    # start empty array for phenotpes
    phens = np.empty([files.shape[0], 62128])


    # Loop for loading and stacking all data
    for i, file in enumerate(files[0]):
        print(file)
        # load the nifty
        img = nib.load(file)
        # Convert to numpy array
        arr = img.get_fdata()
        # Extract upper triangle and flatten it into the the running set of values
        phens[(i-1),:] = arr[np.triu_indices(arr.shape[0])]
    
    
    # convert it to a dataframe
    phens = pd.DataFrame(phens)
    
    # pattern to match for subject ids
    id_pat = "sub-(.*?)_ses"
    # pattern for sesison
    ses_pat = "ses-(.*?)_task"
    # patter for task
    task_pat = "_task-(.*?).pconn"
    
    # seed empty lists
    ids = []
    sess = []
    tasks = []


    #parse id, session, and task
    for i in files[0] : 
        ids.append(re.search(id_pat, i).group(1))
        sess.append(re.search(ses_pat, i).group(1))
        tasks.append(re.search(task_pat, i).group(1))
    
    # Save ID, seession, and task
    phens["IID"] = ids
    phens["ses"] = sess
    phens["task"] = tasks
        
    
    # save if you want
    if (save == True) :
        # add subject ID
        phens["grid"] = files[0].str.split("/", expand = True)[1].str.split("_", expand=True)[0]
        phens.to_csv(save_path, columns = False, index= False, quotes = False)
        
    # return the dataframe
    return(phens)
    
#%%
# test
# os.chdir("/home/christian/Research/Stat_gen/practice_conn_phenos/")
# df = load_extract_niis("/home/christian/Research/Stat_gen/practice_conn_phenos/files")

#%%


