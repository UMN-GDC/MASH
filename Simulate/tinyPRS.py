import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm


N = 1000
M = 250
rng = np.random.default_rng()
Beta0 = np.repeat([5], M)
Xdiff = 0.2
X1 = rng.binomial(2, 0.5 - Xdiff, size = (N, M))
X2 = rng.binomial(2, 0.5 + Xdiff, size = (N, M))
Y1 = np.dot(X1, Beta0) + rng.normal(0, 1, size = N) 
Y2 = np.dot(X2, Beta0) + rng.normal(0, 1, size = N)

pcs = PCA(n_components = 2).fit_transform(np.concatenate([X1, X2]))
sns.scatterplot(x = pcs[:, 0], y = pcs[:, 1], hue = np.concatenate([np.repeat([0], N), np.repeat([1], N)]))
plt.show()


betas = np.zeros(X1.shape[1])
for s in range(X1.shape[1]):
    try :

        betas[s] = sm.regression.linear_model.OLS(Y1,
                                           sm.add_constant(X1[:, s])).fit().params[1]
    except KeyError:
        betas[s] = 0
good =  abs(betas) > 0
PRS = np.dot(X1[:, good], betas[good])
prsModel  = sm.regression.linear_model.OLS(Y1, sm.add_constant(PRS)).fit()
print(prsModel.rsquared)


PRSout = np.dot(X2[:,good], betas[good])
prsOutModel = sm.regression.linear_model.OLS(Y2, sm.add_constant(PRSout)).fit()
print(prsOutModel.rsquared)


