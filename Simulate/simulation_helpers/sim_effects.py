import numpy as np


def sim_effects(rng, X, target_var, frequencyDependent = False, evenlySpaced = False):
    if not evenlySpaced and not frequencyDependent: 
        # simulate effects (beta) from a standard normal distribution making sure that is is the same dimension as the number of columns of X
        beta = rng.standard_normal(X.shape[1]) 
    elif frequencyDependent and not evenlySpaced : 
        # find the average across each row and divide by 2 to get binomial frequency
        freq = np.mean(X, axis = 0)/2
        prop = freq * (1-freq) ** (-1)
        # sample
        beta = rng.normal(np.repeat(0, X.shape[1]), prop, size= X.shape[1])

    elif evenlySpaced and not frequencyDependent : 
        # Make betas evenly spaced around zero for the number of columns X has 
        beta = np.linspace(-1, 1, X.shape[1])
    elif evenlySpaced and frequencyDependent :
        raise ValueError("Cannot have evenly spaced and frequency dependent effects")
    # get the inner product of X and beta
    inner_prod = np.dot(X, beta)
    # rescale inner_product such that the variance is equal target_var
    inner_prod = np.sqrt(target_var / np.var(inner_prod)) * inner_prod

    return inner_prod
