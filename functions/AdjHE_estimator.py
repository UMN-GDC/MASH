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
import statsmodels.api as sm
import timeit
import resource
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from functions.eigenvector_outters import multiple_outer


def AdjHE_estimator(A,data, mp, npc=0, std=False):
    """
    Fucntion for generating heritability estimates from an Adjusted HE standpoint in closed form.

    Parameters
    ----------
    A : numpy array  
        n x n array containing the GRM values.
    data : pandas dataframe
        A dataframe containing all covariates, principal components, and phenotype that has been residualized if necessary.
    mp : int
        1 based index specifying the phenotype to be estimated on.
    npc : int, optional
        number of prinicpal components to adjust for. The default is 0.
    std : bool, optional
        Specify whether or not to standaridize the variables before estimation. The default is False.

    Returns
    -------
    tumple(scalar, scalar)
        h2 - heritability estimate.
        standard error estimate
    """
    # remove identifiers from y for linear algebra 
    y = data[mp]
    # select PC columns 
    PC_cols = [ col.startswith("pc")   for col in data ]
    PCs = data.iloc[:, PC_cols]
    
    trA2 = np.sum(np.multiply(A,A))
    trA = np.sum(np.diag(A))
    n = A.shape[1]

    # If standardized AdjHE is chosen 
    if (std == True) :
        # Standardize the y
        y = (y-np.mean(y))/np.std(y)
    
    yay = np.dot(y.T, np.dot(A, y)).flatten()
    yty = np.dot(y.T, y).flatten()


    # If standardized AdjHE is chosen 
    tn = np.sum(y)**2/n # all 1s PC

    if (npc==0):
        denominator = trA2 - 2*trA + n
        nominator = n - trA + yay - yty
        sigg = n*yay - trA*yty
        sigg = sigg-yay+tn*trA # add 1's
        sige = trA2*yty - trA*yay
        sige = sige-tn*trA2 # add 1's

    else:
        pc = PCs
        s = np.diag(np.dot(pc.T,np.dot(A,pc)))
        b = s - 1
        c = np.dot(y.T, pc)**2 - 1
        nominator = n - trA + yay - yty - np.sum(b*c)
        
        # remove identifiers for linear algebra
        pc = PCs
        pcA = np.dot(pc.T,A)
        pcApc = np.dot(pcA,pc)
        s = np.diag(pcApc) #pciApci
        b = s-1
        t = np.dot(y.transpose(),pc)**2 #ypcipciy
        a11 = trA2 - np.sum(s**2) 
        a12 = trA - np.sum(s)
        b1 = yay - np.sum(s*t)
        b2 = yty - np.sum(t)
        sigg = (n-npc)*b1 - a12*b2
        sigg = sigg.flatten() - yay.flatten() + tn * a12 # add 1's
        sige = a11*b2 - a12*b1
        sige = sige.flatten()-tn*a11 # add 1's
        denominator = trA2 - 2*trA + n - np.sum(b**2)
        
    if std :
        h2= nominator/denominator
        h2 = h2[0]
        
    else:
        h2 = sigg/(sigg+sige)
            
        var_ge = 2/denominator 
        
    return h2,np.sqrt(abs(var_ge))



def create_formula(nnpc, covars, mp):
    """
    Creates a formula for the mean to residualize against. In other words describes the projection space and by extension the residual space.

    Parameters
    ----------
    nnpc : int
        number of PC's included in the projection.
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
    covar_part = " + ".join(covars)
    pc_part = " + ".join(pc_cols)
    
    if (len(covars) == 0) and (len(pc_cols) != 0) :
        RHS = pc_part
    elif (len(covars) != 0) and (len(pc_cols)== 0) :
        RHS = covar_part
    elif (len(covars) == 0) and (len(pc_cols)== 0) :
        RHS = "1"
    else :
        RHS = " + ".join([covar_part, pc_part])


    form = mp + "~ " +  RHS
    # columns
    cols = id_cols + [mp] + covars + pc_cols
    # return the formula and columns
    return(form, cols)
 


def load_n_AdjHE(df, covars, nnpc, mp, ids, GRM_array_nona, std = False):
    """
    Estimates heritability using the efficient AdjHE closed form solution. Takes a dataframe, selects only the
    necessary columns (so that when we do complete cases it doesnt exclude too many samples) residualizes the 
    phenotype, then documents the heritability, standard error and some computer usage metrics.

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
    ids : np array
        np array of subject ids.
    GRM_array_nona : np array
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
    # seed empty result vector
    # result.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]
    # start clock for fitting 
    start_est = timeit.default_timer()
    # create the regression formula and columns for seelcting temporary
    form, cols  = create_formula(nnpc, covars, mp)
    # save a temporary dataframe
    temp = df[cols].dropna()
    # Save residuals of selected phenotype after regressing out PCs and covars
    temp[mp] = smf.ols(formula = form, data = temp, missing = 'drop').fit().resid
    # Potentially could use this to control for random effects
    # smf.mixedlm(formula= form, data = temp, groups=temp["scan_site"])
    # keep portion of GRM without missingess for the phenotypes or covariates
    nonmissing = ids[ids.iid.isin(temp.iid)].index
    GRM_nonmissing = GRM_array_nona[nonmissing,:][:,nonmissing]
    # Get heritability and SE estimates
    h2, se = AdjHE_estimator(A= GRM_nonmissing, data = temp, mp = mp, npc=nnpc, std=std)
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
    print(list(temp.columns))
    print(result["h2"])
    # Return the fit results
    return(pd.DataFrame(result))


#%%


def load_n_MOM(df, covars, nnpc, mp, ids, GRM_array_nona, std = False):
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
    ids : np array
        np array of subject ids.
    GRM_array_nona : np array
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
    # seed empty result vector
    # result.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]
    # start clock for fitting 
    start_est = timeit.default_timer()
    # create the regression formula and columns for seelcting temporary
    form, cols  = create_formula(nnpc, covars, mp)
    # save a temporary dataframe
    temp = df[cols].dropna()
    # Save residuals of selected phenotype after regressing out PCs and covars
    temp[mp] = smf.ols(formula = form, data = temp, missing = 'drop').fit().resid
    # Potentially could use this to control for random effects
    # smf.mixedlm(formula= form, data = temp, groups=temp["scan_site"])
    # keep portion of GRM without missingess for the phenotypes or covariates
    nonmissing = ids[ids.iid.isin(temp.iid)].index
    GRM_nonmissing = GRM_array_nona[nonmissing,:][:,nonmissing]
    # Get heritability and SE estimates
    # Create a matrix of the second order terms
    temp2 = pd.DataFrame({
        # get dependent from outter product of residualized phenotypes
        #"yyt" : np.outer(temp[mp], temp[mp])[np.triu_indices(len(temp[mp]))].flatten(),
        # get GRM
        "A" : GRM_nonmissing[np.triu_indices(len(temp[mp]))].flatten(),
        # get I matrix
        "I" : np.identity(len(temp[mp]))[np.triu_indices(len(temp[mp]))].flatten()})
    # Get outer product of eigen loadings
    temp2 = pd.concat([temp2,
                       multiple_outer(temp, nnpc)], axis = 1)
    temp2[mp] = np.outer(temp[mp], temp[mp])[np.triu_indices(len(temp[mp]))].flatten()
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
    print(list(temp.columns))
    print(result["h2"])
    # Return the fit results
    return(pd.DataFrame(result))



#%%

def load_n_estimate(df, covars, nnpc, mp, ids, GRM_array_nona, std = False, fast=True):
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
    ids : np array
        np array of subject ids.
    GRM_array_nona : np array
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
    
    # Select method of estimation
    if fast == True: 
        print("AdjHE")
        result = load_n_AdjHE(df, covars, nnpc, mp, ids, GRM_array_nona, std = False)

    else: 
        print("OLS")
        result = load_n_MOM(df, covars, nnpc, mp, ids, GRM_array_nona, std = False)
    
    return(pd.DataFrame(result))

#%%
# d = load_n_all_estimate(df, covars, nnpc, "first3000", ids, GRM_array_nona, std = False, fast=False)
