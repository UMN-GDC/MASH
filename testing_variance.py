import numpy as np

# Add the additional alpha +1 so that it can simulate for both alpha = 0 and alapha = -1
homo_eff = rng.normal(np.repeat(0, nSNPs), np.sqrt(sig_g / nSNPs * (2 * prop) ** alpha  * 2 ** (alpha+ 1)), size =  nSNPs)

