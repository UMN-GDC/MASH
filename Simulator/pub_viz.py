import numpy as np
import matplotlib.pyplot as plt


#%% vizualize different alphas on SNP effect size
x = np.linspace(0, 0.25, 20)
#plt.subplots(1,1)
for i in [0,0.25, 0.5, 0.75,1] :
    plt.plot(x, (x * (1-x))**(-i), label = "alpha= " + str(i)),
plt.legend()
plt.xlabel("MAF")
plt.ylabel("SNP effect size")
plt.savefig(fname = "docs/Figures/Sim_features/SNPalphaControl.png", dpi = 300, format = "png")


