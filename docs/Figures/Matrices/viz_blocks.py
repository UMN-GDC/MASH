import numpy as np
from scipy.linalg import block_diag 
import matplotlib.pyplot as plt
ni = int(40)
sites = int(25)



arrays = [np.ones((ni,ni)) for i in range(sites)]
X = block_diag(*arrays)
Y = np.eye(sites* ni)
#%%
fig, ax = plt.subplots(1,2)
ax[0].imshow(X)
ax[1].imshow(Y)
ax[0].set_title("S")
ax[1].set_title("I")
fig.savefig("S_I.png")
