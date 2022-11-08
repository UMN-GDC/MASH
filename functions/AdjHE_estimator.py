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
import timeit
import resource
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from functions.eigenvector_outters import columnwise_outter
from functions.AdjHE_estimator_w_rv import AdjHE_rv_estimator



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




def create_formula(nnpc, covars, mp, RV = None):
    """
    Creates a formula for the mean to residualize against. In other words describes the projection space and by extension the residual space.

    Parameters
    ----------
    nnpc : int
        number of PC's included in model.
    covars : list of integers
        the covariates to be included in the projection.
    mp : int
        phenotype to be projected.

    Returns
    -------
    (form, cols) where form is the formual exceptabule by smf.ols procedure and cols specifies the columns necessary for the projection
    (PCs, covariates, and phenotype).

    """
    # Get indices for ID variables
    id_cols = ["fid", "iid"] 
    # Get the full range of pc columns
    pc_cols = ["pc_" + str(p) for p in range(1, nnpc +1)]
    # Create formula string
    
    if len(covars) != 0:
        RHS = " + ".join(covars)
    else :
        RHS = "1"


    form = mp + "~ " +  RHS
    # columns
    cols = id_cols + [mp] + covars + pc_cols
    if RV != None :
        cols = cols + [RV]

    # return the formula and columns
    return(form, cols)
 


def load_n_AdjHE(df, covars, nnpc, mp, GRM, std = False, RV = None):
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
        result = AdjHE_rv_estimator(A= GRM, df = df, mp = mp,rv=RV, npc=nnpc, std=std)
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
    
    print(list(df.columns))
    print(result["h2"])
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

    print("MOM regression estimator")
    # Get heritability and SE estimates
    # Create a matrix of the second order terms
    temp2 = pd.DataFrame({
        # get dependent from outter product of residualized phenotypes
        #"yyt" : np.outer(temp[mp], temp[mp])[np.triu_indices(len(temp[mp]))].flatten(),
        # get GRM
        "A" : GRM[np.triu_indices(len(df[mp]))].flatten(),
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
    print(list(df.columns))
    print(result["h2"])
    # Return the fit results
    return(pd.DataFrame(result, index = [0]))



#%%

def load_n_estimate(df, covars, nnpc, mp, GRM, std = False, fast=True, RV = None):
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
    ids = df[["fid", "iid"]]
    # seed empty result vector
    # result.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]
    # create the regression formula and columns for seelcting temporary
    form, cols  = create_formula(nnpc, covars, mp, RV)
    # save a temporary dataframe
    temp = df[cols].dropna()
    # Save residuals of selected phenotype after regressing out PCs and covars
    temp[mp] = smf.ols(formula = form, data = temp, missing = 'drop').fit().resid
    # Potentially could use this to control for random effects
    # smf.mixedlm(formula= form, data = temp, groups=temp["scan_site"])
    # keep portion of GRM without missingess for the phenotypes or covariates
    nonmissing = ids[ids.iid.isin(temp.iid)].index
    GRM_nonmissing = GRM[nonmissing,:][:,nonmissing]

    
    # Select method of estimation
    if fast == True: 
        print("AdjHE")
        result = load_n_AdjHE(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)

    else: 
        print("OLS")
        result = load_n_MOM(temp, covars, nnpc, mp, GRM_nonmissing, std = False, RV = RV)
    
    return(pd.DataFrame(result))


    



