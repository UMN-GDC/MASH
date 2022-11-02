#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:02:50 2022

@author: christian
"""

import os
import pandas as pd
import plotly.express as px
from plotly.offline import plot
# os.chdir("/home/christian/Research/Stat_gen/tools/Basu_herit/Genes_n_sims")

#%%

def load_n_plot(file) :
    df = pd.read_csv(file)
    df["h2"] = df.sg / (df.sg + df.ss + df.se)
    df2 = pd.melt(df, id_vars=['sg', 'ss', 'se', 'h2'], value_vars=['Basic_est', 'Site_RE', 'Site_FE'])
    fig = px.violin(df2, x="variable", y="value", color="variable", facet_col="h2", facet_row = "ss")
    fig.update_yaxes(matches=None)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
    fig.update_xaxes(matches=None)
    fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True))
    plot(fig, filename = "Plots/"+ file[0:12] + ".html")

#%%
for f in os.listdir() : 
    if os.path.isfile(f):
        load_n_plot(f)    
        
        
        
        
#%%
df = results
df["h2"] = df.sg / (df.sg + df.ss + df.se)
df["h2"][np.isnan(df.h2)] = 0
df2 = pd.melt(df, id_vars=['sg', 'ss', 'se', 'h2'], value_vars=['Basic_est', 'Site_RE', 'Site_FE'])
fig = px.box(df2, x="variable", y="value", color="variable", facet_col="h2", facet_row = "ss")
fig.update_yaxes(matches=None)
fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
fig.update_xaxes(matches=None)
fig.for_each_xaxis(lambda xaxis: xaxis.update(showticklabels=True))
plot(fig)
#plot(fig, filename = "Plots/"+ file[0:12] + ".html")

