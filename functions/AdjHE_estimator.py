import numpy as np
import statsmodels.api as sm
import timeit
import resource
import numpy as np
import pandas as pd

def AdjHE_estimator(A,data, mp, npc=0, std=False):
    # remove identifiers from y for linear algebra 
    y = data["res" + str(mp)]
    # select PC columns 
    PC_cols = [ col.startswith("PC")   for col in data ]
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

 

def load_n_estimate(data, covar_set, pc_set, pheno_number, fit_number, ids, GRM_array_nona, npc, std):
    result = pd.DataFrame(np.zeros(1, 4))
    result.columns = ["h2", "SE", "Time for analysis(s)", "Memory Usage"]

    # Save temp with just the phenotype we need (I'm sure this can be written given the hints that python returns
    id_cols = (data.columns == "FID") + (data.columns == "IID") 
    pc_cols = data.columns.str.startswith('PC')
    covar_cols = data.columns.str.startswith('Covar')
    pheno_col = (data.columns == 'Pheno_'+ str(pheno_number))
    temp=data.loc[:,id_cols+ pc_cols + covar_cols + pheno_col]
    # drop missing values from both phenos and covariates
    temp = temp.dropna()    
    # Save residuals of selected phenotype after regressing out PCs and covars
    temp["res" + str(pheno_number)] = sm.OLS(endog= temp.loc[:,"Pheno_" + str(pheno_number)], exog= temp.loc[:,temp.columns.str.startswith('PC')+ temp.columns.str.startswith('Covar')]).fit().resid
    # keep portion of GRM without missingess
    nonmissing = ids[ids.IID.isin(temp.IID)].index
    GRM_nonmissing = GRM_array_nona[nonmissing,:][:,nonmissing]
    # resutls from mp pheno
    start_est = timeit.default_timer()
    # Get heritability and SE estimates
    result.iloc[fit_number,0], result.iloc[fit_number,1] = AdjHE_estimator(A= GRM_nonmissing, data = temp, mp = pheno_number, npc=npc, std=std)
    # Get time for each estimate
    result.iloc[fit_number, 2] = timeit.default_timer() - start_est
    # Get memory for each step (This is a little sketchy)
    result.iloc[fit_number, 3] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # return vector of results
    return(result)
   
