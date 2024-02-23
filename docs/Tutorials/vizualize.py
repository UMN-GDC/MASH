import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Estimate.estimators.GCTA_wrapper import GCTA, gcta
from Estimate.estimators.all_estimators import h2Estimation
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
from Simulate.summarizers.genotype_viz import SNPeffectViz, plotClusters
from sklearn.decomposition import PCA

#%% Basic testing for single site and cluster


rng = np.random.default_rng(123)
sim = pheno_simulator(rng = rng, nsubjects= 600)
sim.sim_sites(nsites =1)
sim.sim_pops(nclusts= 3, theta_alleles = [0.1, 0.1], shared = 0.5)
sim.sim_genos(maf = None)
sim.sim_pheno(h2Hom = 0.5, h2Het= [0, 0, 0], alpha = 0, prop_causal =[0,0.5])



# fit and transform 
pca = PCA(n_components=2)
genoPC = pca.fit_transform(sim.genotypes)
clustPC = pca.transform(sim.cluster_frequencies) * 2
ancPC = pca.transform(sim.ancest_freqs.reshape(1,-1))
empCenters = pca.transform(pd.DataFrame(sim.genotypes).assign(cluster = sim.df.subj_ancestries).groupby('cluster').mean())



# plot the genoPCs with the clustPC over
sns.scatterplot(x = genoPC[:,0], y = genoPC[:,1], hue = sim.df.subj_ancestries, palette = 'dark')
sns.scatterplot(x = empCenters[:,0], y = empCenters[:,1], color = "black", marker = "X")
sns.scatterplot(x = clustPC[:,0], y = clustPC[:,1], color = "black")
sns.scatterplot(x = ancPC[:,0], y = ancPC[:,1], color = "black")
plt.text(ancPC[0,0], ancPC[0,1], "Common Ancestry", fontsize = 12)

# label each clustPC
for i in range(sim.nclusts):
    plt.text(clustPC[i,0], clustPC[i,1], f"Ancestry {i}", fontsize = 12)
# Label axes
plt.xlabel('PC1')
plt.ylabel('PC2')

# save the plot
plt.savefig("docs/Figures/samplesNcentersOnPCs.png", dpi = 300)
