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
sim.full_sim(alpha = -1)
sim2 = pheno_simulator(nsubjects= 1000, nSNPs = 1000)
sim2.full_sim(alpha = 0)


MAF = 0.5 - abs(0.5 - sim.genotypes[:, sim.causals].mean(axis= 0)/2)
EFF =abs(sim.SNP_effects)
EFF = EFF/max(EFF)
EFF2 = abs(sim2.SNP_effects)
EFF2=  EFF2/max(EFF2) * 0.5
MAF2 = 0.5 - abs(0.5 - sim2.genotypes[:, sim2.causals].mean(axis= 0)/2)
plt.scatter(x = MAF, y=  EFF, label = "alpha = -1",alpha = 0.5)
plt.scatter(x = MAF2, y=  EFF2, label = "alpha = 0", alpha = 0.5)
plt.legend() 
plt.xlabel("MAF")
plt.ylabel("SNP effect")
plt.savefig(fname = "docs/Papers/Toolbox_paper/figures/simmedSNPeffects.png", dpi = 300, format = "png")


