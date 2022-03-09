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

def load_extract_niis(files, IDs, save = False, save_path = None) :
    # load the list of file paths
    files = pd.read_csv(files, header=None)

    # start empty array for phenotpes
    phens = np.empty([files.shape[0], 62128])


    # Loop for loading and stacking all data
    for i, file in enumerate(files[0]):
        # print(file)
        # load the nifty
        img = nib.load(file)
        # Convert to numpy array
        arr = img.get_fdata()
        # Extract upper triangle and flatten it into the the running set of values
        phens[(i-1),:] = arr[np.triu_indices(arr.shape[0])]
    
    
    # convert it to a dataframe
    phens = pd.DataFrame(phens)
    
    # pattern to match for subject ids
    id_pat = "sub-(.*?)/ses"
    # pattern for sesison
    ses_pat = "ses-(.*?)_task"
    # patter for task
    task_pat = "_task-(.*?).pconn"
    
    # seed empty lists
    iids = []
    sess = []
    tasks = []


    #parse id, session, and task
    for i in files[0] :
        # select IID from the file name
        iid =  re.search(id_pat, i).group(1)
        # Add a - to make names consistent across data types
        iid = iid[:4] + "_" + iid[4:]
        # Add to cumulative iid column
        iids.append(iid)
        sess.append(re.search(ses_pat, i).group(1))
        tasks.append(re.search(task_pat, i).group(1))
    
    # Save ID, seession, and task
    phens.insert(loc=0, column='IID', value=iids)    
    #phens["ses"] = sess
    #phens["task"] = tasks
    phens = pd.merge(IDs, phens, on="IID") 
    # Save ID's and FID's to reenter in front of dataframe
    IIDs = phens.IID
    FIDs = phens.FID
    phens.drop(["FID", "IID"], axis = 1)
    
    #phens.insert(loc=0, column= "IID", value= IIDs)
    #phens.insert(loc=0, column= "FID", value= FIDs)

    
    # save if you want
    if (save == True) :
        phens.to_csv(save_path, columns = False, index= False, quotes = False)
        
    # return the dataframe
    return(phens)
    
#%%
# test
# os.chdir("/home/christian/Research/Stat_gen/practice_conn_phenos/")
# df = load_extract_niis("/home/christian/Research/Stat_gen/practice_conn_phenos/files")

#%%


