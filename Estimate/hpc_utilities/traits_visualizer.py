#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:03:07 2022
Christian Coffman: coffm049@umn.edu
Created 2022-08-22
Last Updated 2022-08-22
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

#%% Create a small test dataframe for trouble shooting

# df = pd.DataFrame({"d1" : ["a", "b", "c"] * 10,
#                     "d2" : np.repeat(["a", "b", "c"], 10),
#                     "d3" : np.repeat("a", 30),
#                     "c1" : np.random.rand(30,),
#                     "c2" : np.random.rand(30,)})


#%%
def disc_plot(df, d1, d2, out) :
    """
    Plots the percentages with which discrete variable d2 is observed across discrete variable d1.

    Parameters
    ----------
    df : DataFrame
        Long format dataframe contianing a row for each observation of the dataset
    d1 : string
        String specifying the column name for the discrete variable to be depicted across the x-axis.
    d2 : string
        String specifying the column name for the discrete variable to count per each level of d1 and comprise the different colors on the bargraph.
    out : string
        specify folder to save image in
    
    Returns
    -------
    None. Saves image to specified output folder


    """
    
    # Get counts of data by the specified discrete variables
    
    test = (
            df.value_counts([d1, d2])
            .rename("var")
            )
    
    # Get the totals by the first specified discrete variable
    groupsums = test.groupby(d1).sum("var")

    # Standardize each group by the total within that group
    test = test / groupsums * 100

    # plot it
    test = (
            test.reset_index()[[d1, d2, "var"]]
            .pivot(index = d1, columns = d2)
            )


    # plot it
    p = test.plot(kind = "bar", stacked = True, y="var")
    p.set(xlabel = d1, ylabel = "Percentage")
    
    # save plot
    fig = p.get_figure()
    fig.savefig(out + d2 + "_vs_" + d1)
    plt.close()


def cont_plot(df, d1, c, out) :
    """
    boxplots for the continuous trait across the discrete trait in question.

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe containing one row per observation.
    d1 : string
        specify the trait to use along the x-axis for the boxplot.
    c : string
        specify the trait to compare across all levels of d1.
    out : string
        specify folder to save image in
    
    Returns
    -------
    None. Saves image to specified output folder

    """
    p = sns.boxplot(data= df, x = d1, y = c)
    #p.set(xlabel= d1, ylabel = c)
    
    fig = p.get_figure()
    fig.savefig(out + c + "_vs_" + d1)
    plt.close()
    

def covs_vs_cov_of_interest(df, xvar, other_vars, out) :
    """
    Plot relationship of variables to variable of interest (xvar), boxplots for continouosu variables, and stacked bargrpahs for discrete.
    This will be especially valuable once 

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe containing all necessary variables with one observation per row (long format).
    xvar : string
        Variable of interest to compare to all other variables.
    other_vars: list of integers
        Specifies which covariates are of interest
    out : string
        specifies file path for analysis directory.

    Returns
    -------
    None.

    """
    # Make a directory for var relations if it doesn't already exist
    os.makedirs(out + "/Results/Var_relations", exist_ok=True)
    
    for col in df.columns:
    #for col in df.columns[[np.array(other_vars) + 2]]:
        # test if column is the random variable the user specified
        if col == xvar : (
                # This is the random variable
                )
        # test if column is FID or IID
        elif col in ["fid", "iid"] : (
                # These are unique identifiers
                )
        else :
                # test if data has more than 30 unique values. If so consider it continuous, 
            if len(set(df[col])) >= 30 :
                cont_plot(df, xvar, col, out + "/Results/Var_relations/")
            # else treat it as discrete
            else :
                disc_plot(df, xvar, col, out + "/Results/Var_relations/")
