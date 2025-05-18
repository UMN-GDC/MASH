import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from Simulate.simulation_helpers.Sim_generator import pheno_simulator 

def reg_coef(x,y,label=None,color=None,**kwargs):
    ax = plt.gca()
    r,p = pearsonr(x,y)
    ax.annotate('r = {:.2f}'.format(r), xy=(0.5,0.5), xycoords='axes fraction', ha='center')
    ax.set_axis_off()

#%%
Similar = False 
sim = pheno_simulator()
sim.sim_sites(nsites = 1)
sim.sim_pops(nclusts = 2 )
sim.sim_genos()

titleString = "SNPeffects_"
if Similar : 
    sim.sim_pheno(h2Hom = 0.8, h2Het = [0.05,0.05])
    titleString += "Similar.png"
else :
    sim.sim_pheno(h2Hom = 0.05, h2Het = [0.8,0.8])
    titleString += "Different.png"


het2 = sim.het_eff

for i in range(sim.het_eff.shape[1]) : 
    het2[:,i] += sim.homo_eff


g = sns.PairGrid(pd.DataFrame(het2))
g.map_diag(sns.distplot)
g.map_lower(sns.regplot)
g.map_upper(reg_coef)
g.savefig(titleString)

