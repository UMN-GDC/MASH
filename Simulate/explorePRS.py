import sys
import os
import logging
import subprocess
import numpy as np
import pandas as pd 
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns

#%%
rng = np.random.default_rng(123)
sim = pheno_simulator(rng = rng, nsubjects= 100)
sim.sim_sites(nsites =1)
sim.sim_pops(nclusts= 2, theta_alleles=1e-6)
sim.sim_genos()
X = StandardScaler().fit_transform(sim.genotypes)
pcs = PCA(n_components = 2).fit_transform(sim.genotypes)
sns.scatterplot(x = pcs[:, 0], y = pcs[:, 1], hue = sim.df["subj_ancestries"])
plt.show()

sim.sim_pheno(h2Hom = 0.75, h2Het= [0, 0], alpha = 0)
df1= sim.df.query("subj_ancestries==0")
scaler = StandardScaler()
G1 = scaler.fit_transform(sim.genotypes[df1.index, :])
df2 = sim.df.query("subj_ancestries==1")
G2 = scaler.fit_transform(sim.genotypes[df2.index, :])

#%% PRS process
betas = np.zeros(G1.shape[1])
for s in range(G1.shape[1]):
    try :

        betas[s] = sm.regression.linear_model.OLS(df1["Y0"],
                                           sm.add_constant(G1[:, s])).fit().params["x1"]
    except KeyError:
        betas[s] = 0
good =  abs(betas) > 0
PRS = np.dot(G1[:, good], betas[good])
prsModel  = sm.regression.linear_model.OLS(df1["Y0"], sm.add_constant(PRS)).fit()
print(prsModel.rsquared)


PRSout = np.dot(G2[:,good], betas[good])
prsOutModel  = sm.regression.linear_model.OLS(df2["Y0"], sm.add_constant(PRSout)).fit()
print(prsOutModel.rsquared)

plt.plot(PRS, df1["Y0"], "o")
plt.plot(PRSout, df2["Y0"], "o")
plt.show()
