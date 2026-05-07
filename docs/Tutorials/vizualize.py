import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from Simulate.simulation_helpers.Sim_generator import pheno_simulator

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

print("PCA completed successfully")
print(f"Genotype PCs shape: {genoPC.shape}")
print(f"Cluster centers shape: {clustPC.shape}")
print(f"Ancestry center shape: {ancPC.shape}")