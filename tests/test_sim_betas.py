import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
rng = np.random.default_rng(123)

#%%
sim = pheno_simulator(rng = rng)
sim.sim_sites()
sim.sim_pops(nclusts= 2)
sim.sim_genos()

# Similar or sissimilar SNP effects 
Similar = False 
if Similar :
    h2Hom = 0.5
    h2Het = [0.01,0.01]
else : 
    h2Hom = 0.01
    h2Het = [0.5,0.5]

sim.sim_pheno(h2Hom = h2Hom, h2Het = h2Het)
eff1 = sim.homo_eff + sim.het_eff[:,0]
eff2 = sim.homo_eff + sim.het_eff[:,1]
#%% plot eff1 as a scatter plot with index as x axis 
fig, axes = plt.subplots(2, 1)
fig.suptitle(f"SNP effects across populations (h2Hom= {str(h2Hom)}, h2Het = {str(h2Het)})")
sns.scatterplot(ax = axes[0], x = range(len(eff1)), y=eff1, c = "blue")
sns.scatterplot(ax = axes[0], x = range(len(eff2)), y= eff2, c = "red")
sns.lineplot(ax = axes[0], x = range(len(eff1)), y = abs(eff1-eff2), c= "black")
sns.kdeplot(abs(eff1 - eff2), ax = axes[1])
plt.savefig(f"docs/Figures/Phenotypes/SNPEffectsSimilar_{str(Similar)}.png")
del sim
