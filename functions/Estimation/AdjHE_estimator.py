#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AdjHE estimator guts
Christian Coffman: coffm049@umn.edu
Created 2022-05-26
Last Updated 2022-05-26
"""

##############################################################
# These helpers stores all of the functions necessary for the 
# guts of the AdjHE estimator
##############################################################

import numpy as np
from numpy import matrix
from scipy.linalg import block_diag
import timeit
import resource
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from functions.Estimation.Estimate_helpers import create_formula, columnwise_outter



def AdjHE_estimator(A, df, mp, npc=0, std=False):
    """
    Fucntion for generating heritability estimates from an Adjusted HE standpoint in closed form.
    Parameters
    ----------
    A : numpy array  
        n x n array containing the GRM values.
    df : pandas dataframe
        A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
    mp : int
        1 based index specifying the phenotype to be estimated on.
    npc : int, optional
        number of prinicpal components to adjust for. The default is 0.
    std : bool, optional
        Specify whether or not to standaridize the variables before estimation. The default is False.
    Returns
    -------
    tumple(ndarray, scalar)
        variance parameter estimates.
        standard error estimate
    """
    y = np.array(df[mp])
    A = np.array(A)

    if npc ==0 :
        Sjs = 0
        Tjs =0
    else :
        # Grab the PCs
        PC_cols = [ col.startswith("pc")   for col in df ]
        # Select specified number of Pcs
        PC = np.array(df.iloc[:,PC_cols])[:, :npc]
        # Stack the PCs for easier matrix multiplications
        PCs = np.reshape(PC.T, newshape=(PC.shape[1], PC.shape[0], 1), order = "C")
        # transpose the PCs
        PCsT = np.reshape(PC.T, newshape=(PC.shape[1], 1, PC.shape[0]), order = "C")
    
    
        # Calculate outer products of eigenvectors
        PPt = np.matmul(PCs, PCsT)
    
        # Calculate tj's
        Tjs = np.matmul(np.matmul(y, PPt), y.T)
    
        # Calculated Sj's
        Sjs = np.matmul(np.matmul(PCsT, A), PCs).flatten()

    # Compute elements of regression matrix
    n = A.shape[0]
    
    # Could square A faster with following A^2 = VD^2 V', but need to pass in eigenvalues D
    trA2 = np.trace(np.linalg.matrix_power(A, 2))
    trA = np.trace(A)
    
    topleft = trA2 - np.sum(Sjs ** 2)
    offdiag= trA - np.sum(Sjs)
    bottomright = n- npc
        
    # Solve the regression problem
    XXinv = np.linalg.inv(np.matrix([[topleft, offdiag],
                                     [offdiag, bottomright]]))
    yay = np.matmul(y.T, np.matmul(A, y))
    yty = np.inner(y,y)
    Ycol= np.array([yay - np.sum(Tjs * Sjs),
                    yty - np.sum(Tjs)])
    sigmas = np.array(np.matmul(XXinv, Ycol)).flatten()
    
    
    # Get the variance
    denominator = trA2 - 2*trA + n - np.sum(Sjs**2)
    var_h2 = 2/ denominator
    results = {"sg" : sigmas[0], "ss": 0, "se" : sigmas[1], "var(sg)" : var_h2}
    return results

def AdjHE_rv_estimator(A,df, mp, rv, npc=0, std=False) :
    """
    Estimate the heritability of the presence of an additional random effect.

    Parameters
    ----------
    A : numpy array  
        n x n array containing the GRM values.
    df : pandas dataframe
        A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
    mp : int
        1 based index specifying the phenotype to be estimated on.
    rv : string
        specifying the name of the column to be used as a random variable.
    npc : int, optional
        number of prinicpal components to adjust for. The default is 0.
    std : bool, optional
        Specify whether or not to standaridize the variables before estimation. The default is False.

    Returns
    -------
    tuple(scalar, scalar)
        h2 - heritability estimate.
        standard error estimate
    """
    # Reorder df by the random variable
    # then reorder the GRM to match
    df = df.reset_index().drop("index", axis = 1)
    df = df.sort_values(rv).dropna(subset= [rv])
    A = np.matrix(A[df.index,:][:,df.index])
    n = A.shape[0]
    
    df["Intercept"] = 1
    X_sites= df[rv]

    # Get dummies for categoricals if they exist
    X = np.matrix(pd.get_dummies(df.drop(["FID", "IID", rv], axis = 1),  drop_first = True))
    y = np.matrix(df[mp])

    # Create S similarity matrix 
    site, sizes= np.unique(df[rv], return_counts = True)
    # Construct the block diagonal
    diags = [np.ones((size,size)) for size in sizes]
    S = np.matrix(block_diag(*diags))
    # Standardize S
    # S = (S - S.mean(axis = 1))/ S.std(axis = 1)

    # diags = [np.ones((size,size))* size for size in sizes]
    # S2 = np.matrix(block_diag(*diags) )
    
    # Construct the orthogonal projection matrix Q utilizing QR decomposition
    q, r = np.linalg.qr(X)
    Q = np.identity(n) - X.dot(np.linalg.inv(r).dot(q.T))
    Q = np.matrix(Q)
        
    # Compute elements of 3x3 matrix
    QSQ = Q * S * Q
    
    # find necessary traces
    trA2 = np.trace(A ** 2)
    trQSQA = np.trace(QSQ * A)
    trA = np.trace(A)
    trQSQ = np.trace(QSQ)
    trQSQQSQ = np.trace(QSQ * QSQ)

    # Find inverse 
    XtXm1 = np.linalg.inv(np.matrix([[trA2, trQSQA, trA],
                 [trQSQA, trQSQQSQ, trQSQ],
                 [trA, trQSQ, n]]))

    youter = np.matrix(np.outer(y,y))
    trAY = np.trace(A * youter)
    trQSQY = np.trace(QSQ * youter)
    trYout =  np.trace(youter)
    # Possible that the y's will need to account for prinicpal componetns in future real data cases
    sigmas = XtXm1.dot(np.matrix([[trAY], [trQSQY], [trYout]]))
    
    results = {"sg" : sigmas[0,0], "ss": sigmas[1,0], "se" : sigmas[2,0], "var(sg)" : 0}
    
    # return heritability estimate
    return results



 


def load_n_AdjHE(df, covars, nnpc, mp, GRM, std = False, RV = None, homo = True):
    """
    Estimates heritability using the efficient AdjHE closed form solution. Takes a dataframe, selects only the
    necessary columns (so that when we do complete cases it doesnt exclude too many samples) residualizes the 
    phenotype, then documents the heritability, standard error and some computer usage metrics.

    Parameters
    ----------
    df : pandas dataframe
        dataframe contianing phenotype, covariates, an prinicpal components.
    covars : list of strings
        list of variable names to include in the resiudalization.
    nnpc : int
        number of pcs to include.
    mp : string
        which phenotype to estiamte on.
    GRM : np array
        the GRM with missingness removed.
    std : bool, optional
        specifying whether standarization happens before heritability estimation. The default is False.
    RV : string, optional
        Varible to control for as a random effect, if applicable


    Returns
    -------
    pandas dataframe containing:
        - heritability estimate
        - standard error the estimate
        - the phenotype
        - the number of pcs included
        - The covarites included 
        - time for analysis
        - maximum memory usage
    """
    # start clock for fitting 
    start_est = timeit.default_timer()
    
    # Get heritability and SE estimates from appropriate estimator
    if RV == None :
        result = AdjHE_estimator(A= GRM, df = df, mp = mp, npc=nnpc, std=std)
    else :
        if homo == True :
            result = AdjHE_rv_estimator(A= GRM, df = df, mp = mp,rv=RV, npc=nnpc, std=std)
        else : 
            result = AdjHE_rv_estimator_homo(A= GRM, df = df, mp = mp,rv=RV, npc=nnpc, std=std)
    # Get time for each estimate
    t = timeit.default_timer() - start_est
    # Get memory for each step (in Mb) (This is a little sketchy)
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000
    # Save the formula for the control variables and results
    result = {"h2" : result["sg"] / (result["sg"] + result["ss"] + result["se"]),
              "SE" : np.sqrt(result["var(sg)"]),
              "Pheno" : mp,
              "PCs" : nnpc,
              "Covariates" : "+".join(covars),
              "Time for analysis(s)" : t,
              "Memory Usage" : mem}
    if result["h2"] < 0 :
        result["h2"] = 0
    elif result["h2"] >1 :
        result["h2"] = 1
    
    # Return the fit results
    return(pd.DataFrame(result, index = [0]))


#%%


def load_n_MOM(df, covars, nnpc, mp, GRM, std = False, RV  = None):
    """
    Estimates heritability, but solves a full OLS problem making it slower than the closed form solution. Takes 
    a dataframe, selects only the necessary columns (so that when we do complete cases it doesnt exclude too many samples)
    residualizes the phenotype, then documents the heritability, standard error and some computer usage metrics.

    Parameters
    ----------
    df : pandas dataframe
        dataframe contianing phenotype, covariates, an prinicpal components.
    covars : list of int
        list of integers specifying which covariates to include in the resiudalization.
    nnpc : int
        number of pcs to include.
    mp : int
        which phenotype to estiamte on.
    GRM : np array
        the GRM with missingness removed.
    std : bool, optional
        specifying whether standarization happens before heritability estimation. The default is False.

    Returns
    -------
    pandas dataframe containing:
        - heritability estimate
        - standard error the estimate
        - the phenotype
        - the number of pcs included
        - The covarites included 
        - time for analysis
        - maximum memory usage
    """
    # start clock for fitting 
    start_est = timeit.default_timer()

    # Get heritability and SE estimates
    # Create a matrix of the second order terms
    temp2 = pd.DataFrame({
        # get dependent from outter product of residualized phenotypes
        #"yyt" : np.outer(temp[mp], temp[mp])[np.triu_indices(len(temp[mp]))].flatten(),
        # get GRM
        "A" : np.array(GRM)[np.triu_indices(len(df[mp]))],
        # get I matrix
        "I" : np.identity(len(df[mp]))[np.triu_indices(len(df[mp]))].flatten()})
    # Get outer product of eigen loadings
    temp2 = pd.concat([temp2,
                       columnwise_outter(df, nnpc)], axis = 1)
    temp2[mp] = np.outer(df[mp], df[mp])[np.triu_indices(len(df[mp]))].flatten()
    # Fit the model
    model = sm.OLS(endog = temp2[mp], exog = temp2.drop(mp, axis = 1)).fit()
    sg = model.params["A"]
    stot = model.params.sum()
    h2 = sg / stot
    se = (sg**2/ stot**2) * ((model.bse["A"]**2 / sg**2) + (model.bse ** 2 / stot**2))
    # botom_var = (model.bse**2).sum()

    # Get time for each estimate
    t = timeit.default_timer() - start_est
    # Get memory for each step (in Mb) (This is a little sketchy)
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000
    # Save the formula for the control variables and results
    result = {"h2" : h2,
              "SE" : se,
              "Pheno" : mp,
              "PCs" : nnpc,
              "Covariates" : "+".join(covars),
              "Time for analysis(s)" : t,
              "Memory Usage" : mem}
    # Return the fit results
    return(pd.DataFrame(result, index = [0]))


# def AdjHE_rv_estimator_homo(A,df, mp, rv, npc=0, std=False) :
#     """
#     Estimate the heritability of the presence of an additional random effect.

#     Parameters
#     ----------
#     A : numpy array  
#         n x n array containing the GRM values.
#     df : pandas dataframe
#         A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
#     mp : int
#         1 based index specifying the phenotype to be estimated on.
#     rv : string
#         specifying the name of the column to be used as a random variable.
#     npc : int, optional
#         number of prinicpal components to adjust for. The default is 0.
#     std : bool, optional
#         Specify whether or not to standaridize the variables before estimation. The default is False.

#     Returns
#     -------
#     tuple(scalar, scalar)
#         h2 - heritability estimate.
#         standard error estimate
#     """
#     # Reorder df by the random variable
#     # then reorder the GRM to match
#     df = df.reset_index().drop("index", axis = 1)
#     df = df.sort_values(rv).dropna(subset= [rv])
#     A = np.matrix(A[df.index,:][:,df.index])
#     n = A.shape[0]
    
#     df["Intercept"] = 1
#     X_sites= df[rv]

#     # Get dummies for categoricals if they exist
#     X = np.matrix(pd.get_dummies(df.drop(["FID", "IID", rv], axis = 1),  drop_first = True))
#     y = np.matrix(df[mp])

#     # Create S similarity matrix 
#     site, sizes= np.unique(df[rv], return_counts = True)
#     # Construct the block diagonal
#     diags = [np.ones((size,size)) for size in sizes]
#     S = np.matrix(block_diag(*diags))
#     # Standardize S
#     # S = (S - S.mean(axis = 1))/ S.std(axis = 1)

#     # diags = [np.ones((size,size))* size for size in sizes]
#     # S2 = np.matrix(block_diag(*diags) )
    
#     # Construct the orthogonal projection matrix Q utilizing QR decomposition
#     q, r = np.linalg.qr(X)
#     Q = np.identity(n) - X.dot(np.linalg.inv(r).dot(q.T))
#     Q = np.matrix(Q)
        
#     # Compute elements of 3x3 matrix
#     QSQ = Q * S * Q
    
#     # find necessary traces
#     trA2 = np.trace(A ** 2)
#     trQSQA = np.trace(QSQ * A)
#     trQSQQSQ = np.trace(QSQ * QSQ)

#     youter = np.matrix(np.outer(y,y))

#     denom = trA2 * trQSQQSQ - trQSQA **2 
#     top = (A * trQSQQSQ - QSQ * trQSQA) * youter / denom
#     bottom = (QSQ* trA2 - A * trQSQA) * youter / denom

    
#     results = {"sg" : top, "ss": bottom, "se" : 0, "var(sg)" : 0}
    
#     # return heritability estimate
#     return results
