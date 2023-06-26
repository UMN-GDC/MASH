import numpy as np
import matplotlib.pyplot as plt
from Simulator.simulation_helpers.Sim_generator import pheno_simulator

#%% vizualize different alphas on SNP effect size
x = np.linspace(0, 0.25, 20)
#plt.subplots(1,1)
for i in [0,0.25, 0.5, 0.75,1] :
    plt.plot(x, (x * (1-x))**(-i), label = "alpha= " + str(i)),
plt.legend()
plt.xlabel("MAF")
plt.ylabel("Variance of SNP effect")
plt.savefig(fname = "docs/Papers/Toolbox_paper/figures/SNPalphaControl.png", dpi = 300, format = "png")
plt.clf()

#%% vizualize example simulated effects

sim = pheno_simulator(nsubjects= 1000, nSNPs = 1000)
sim.full_sim()

MAF = 0.5 - abs(0.5 - sim.genotypes[:, sim.causals].mean(axis= 0)/2)
plt.scatter(x = MAF, y=  abs(sim.SNP_effects))
plt.xlabel("MAF")
plt.ylabel("SNP effect")
plt.savefig(fname = "docs/Papers/Toolbox_paper/figures/simmedSNPeffects.png", dpi = 300, format = "png")

