import numpy as np
import statsmodels.api as sm
import timeit
import resource
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

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

def create_formula(nnpc, covars, mp):
    # Get indices for ID variables
    id_cols = ["FID", "IID"] 
    # Get the full range of pc columns
    pc_cols = ["PC_" + str(p) for p in range(1, nnpc +1)]
    # And pheno string
    pheno_col ='Pheno_'+ str(mp)
    # Create formula string
    form = pheno_col + " ~ " + " + ".join(covars) + " + " +  " + ".join(pc_cols)
    # columns
    cols = id_cols + [pheno_col] + covars + pc_cols
    # return the formula and columns
    return(form, cols)
 


def load_n_estimate(df, covars, nnpc, mp, ids, GRM_array_nona, std = False):
    # seed empty result vector
    result = pd.DataFrame(np.zeros((1, 7)))
    result.columns = ["h2", "SE", "Pheno", "PCs", "Time for analysis(s)", "Memory Usage", "formula"]
    # start clock for fitting 
    start_est = timeit.default_timer()
    # create the regression formula and columns for seelcting temporary
    form, cols  = create_formula(nnpc, covars, mp)
    # save a temporary dataframe
    temp = df[cols].dropna()
    # Save residuals of selected phenotype after regressing out PCs and covars
    temp["res" + str(mp)] = smf.ols(formula = form, data = temp, missing = 'drop').fit().resid
    # Potentially could use this to control for random effects
    # smf.mixedlm(formula= form, data = temp, groups=temp["scan_site"])
    # mod = smf.ols(formula='Pheno_' + str(mp) '~ Literacy + Wealth + Region', data=df)
    # keep portion of GRM without missingess
    nonmissing = ids[ids.IID.isin(temp.IID)].index
    GRM_nonmissing = GRM_array_nona[nonmissing,:][:,nonmissing]
    # Get heritability and SE estimates
    result.iloc[0,0], result.iloc[0,1] = AdjHE_estimator(A= GRM_nonmissing, data = temp, mp = mp, npc=nnpc, std=std)
    result.iloc[0, 2] = mp
    result.iloc[0, 3] = nnpc
    # Get time for each estimate
    result.iloc[0, 4] = timeit.default_timer() - start_est
    # Get memory for each step (in Mb) (This is a little sketchy)
    result.iloc[0, 5] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000
    # Save the formula for the control variables
    result.iloc[0, 6] = form
    print(temp.columns)
    print(result.iloc[0,0])
    # Return the fit results
    return(result)




   
