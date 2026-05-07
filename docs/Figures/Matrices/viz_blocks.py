import numpy as np
from scipy.linalg import block_diag 
ni = int(40)
sites = int(25)

arrays = [np.ones((ni,ni)) for i in range(sites)]
X = block_diag(*arrays)
Y = np.eye(sites* ni)

print(f"Matrix X shape: {X.shape}")
print(f"Matrix Y shape: {Y.shape}")