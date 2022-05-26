#! /usr/bin/env python3

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
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

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
    scalar, scalar
        h2 - heritability estimate.
        standard error estimate

    """
    # remove identifiers from y for linear algebra 
    y = data[mp]
    # select PC columns 
    PC_cols = [ col.startswith("pc")   for col in data ]
    PCs = data.iloc[:, PC_cols]
    # If standardized AdjHE is chosen 
    if (std == True) :
        # Standardize the y
        std_y = (y-np.mean(y))/np.std(y)
        
        
        trA = np.sum(np.diag(A))
        trA2 = np.sum(np.multiply(A,A))
        n = A.shape[1]
        yay = np.dot(std_y.T, np.dot(A,std_y)).flatten()
        yty = np.dot(std_y.T, std_y).flatten()
        if (npc==0):
            denominator = trA2 - 2*trA + n
            nominator = n - trA + yay - yty
        else:
            pc = PCs
            s = np.diag(np.dot(pc.T,np.dot(A,pc)))
            b = s - 1
            c = np.dot(std_y.T, pc)**2 - 1
            denominator = trA2 - 2*trA + n - np.sum(b**2)
            nominator = n - trA + yay - yty - np.sum(b*c)
        h2 = nominator/denominator
        h2 = h2[0]
        var_ge = 2/denominator
        #    tau = n/nmarkers
        #    b1 = (1-np.sqrt(tau))**2
        #    b2 = (1+np.sqrt(tau))**2
        #    r = b2-b1
        #    a1 = h2-1
        #    a2 = 1-2*h2
        #    trace_A2_MP = 0.5*(r+2*b1)*n
        #    trace_A3_MP = (5/16*r**2+b1*b2)*n
        #    trace_A4_MP = (7*r**3+30*b1*r**2+48*b1**2*r+32*b1**3)/32*n
        #    if (npc==0):
        #    #    var_MP = 2/denominator
        #        var_ge = 2/denominator
        #    else:
        #        trace_A_MP = trA - np.sum(s)
        #        a = denominator
        #    #    var_MP=2/a**2*(h2**2*trace_A4_MP+(n-npc)*a1**2+(a2**2+2*h2*a1)*trace_A2_MP+2*a1*a2*trace_A_MP+2*h2*a2*trace_A3_MP)
        #        var_ge = 2/a

        
    else :
        
    # else we solve the unstandardized version
        trA2 = np.sum(np.multiply(A,A))
        trA = np.sum(np.diag(A))

        n = A.shape[1]
        yay = np.dot(y.T, np.dot(A,y)).flatten()
        yty = np.dot(y.T, y).flatten()
        tn = np.sum(y)**2/n # all 1s PC
        if (npc==0):
            sigg = n*yay - trA*yty
            sigg = sigg-yay+tn*trA # add 1's
            sige = trA2*yty - trA*yay
            sige = sige-tn*trA2 # add 1's
            denominator = trA2 - 2*trA + n
        else:
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
        h2 = sigg/(sigg+sige)
        var_ge = 2/denominator
    return h2,np.sqrt(var_ge)

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
    form = mp + "~ " + " + ".join(covars) + " + " +  " + ".join(pc_cols)
    # columns
    cols = id_cols + [mp] + covars + pc_cols
    # return the formula and columns
    return(form, cols)
 


def load_n_estimate(df, covars, nnpc, mp, ids, GRM_array_nona, std = False):
    """
    Takes a dataframe, selects only the necessary columns (so that when we do complete cases it doesnt exclude too many samples)
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




   
