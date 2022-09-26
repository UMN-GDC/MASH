import pandas as pd
import numpy as np
from struct import unpack, calcsize
import timeit
from functions.loading_extracting_niis import load_extract_niis

#%%



# The following code was adapted from the GCTA creators website for loading binary GRM's
# https://yanglab.westlake.edu.cn/software/gcta/#MakingaGRM

def sum_n_vec(n):
    s = [int(0)] * n
    for i in range(n):
        s[i] = int(((i + 1) * (i + 2) / 2) - 1)
    return s


def ReadGRMBin(prefix, AllN = False):
    print("Reading GRM: ", prefix)

    # Time reading the GRM and other data
    start_read = timeit.default_timer()

    BinFileName  = prefix + ".grm.bin"
    NFileName = prefix + ".grm.N.bin"
    IDFileName = prefix + ".grm.id"
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
    entry_format = 'f' # N is stored as a float in the binary file
    entry_size = calcsize(entry_format)
    ## Read IDs
    ids = pd.read_csv(IDFileName, sep = '\t', header = None)
    ids_vec = ids.iloc[:,1]
    n = len(ids.index)
    ids_diag = ['NA' for x in range(n)]
    n_off = int(n * (n - 1) / 2)
    ## Read relatedness values
    grm = np.fromfile(BinFileName, dtype = dt)
    ## Read number of markers values
    if AllN:
        N = np.fromfile(NFileName, dtype = dt)
    else:
        with open(NFileName, mode='rb') as f:
            record = f.read(entry_size)
            N = unpack(entry_format, record)[0]
            N = int(N)
    i = sum_n_vec(n)
    ids = ids.rename(columns={0: "fid", 1: "iid"})
    ids = ids.dropna()
    # ids["FID"] = ids.FID.astype(int)
    n_phen_nona = ids.shape[0]
    n_phen_nona = grm[i].size
    GRM_array_nona = np.zeros((n_phen_nona, n_phen_nona))
    val = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids,'N':N, "n_phen_nona" : n_phen_nona}
    return val

def build_grm(G) :
    # Get specific detials about the GRM
    ids = G['id']
    k = G['n_phen_nona']
    GRM = np.zeros((k , k ))
    GRM[np.diag_indices(k)] = G['diag']
    ############################### reconstruct GRM
    # 
    temp_i = 0
    temp = 0
    # k= args.k

    l = list(range(k, k, k))
    l.append(k)
    for i in l:
        cor = multirange(range(temp_i, i))
        GRM[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
        GRM.T[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
        temp = temp + len(cor['b'])
        del(cor)
        temp_i = i
    ################################
    return GRM, ids




# This allows us to read in multiple partial GRM's into one full GRM
def multirange(counts):
    counts = np.asarray(counts)
    # Remove the following line if counts is always strictly positive.
    counts = counts[counts != 0]
    counts1 = counts[:-1]
    reset_index = np.cumsum(counts1)
    incr = np.ones(counts.sum(), dtype=int)
    incr[0] = 0
    incr[reset_index] = 1 - counts1
    # Reuse the incr array for the final result.
    incr.cumsum(out=incr)
    val = {'a':incr,'b':np.repeat(counts,counts)}
    return val




# Read covariates, PC's, and phenotype all at once
def load_data(pheno_file=None, cov_file=None, PC_file= None) :
    
    # Need to specify at least one of pheno or cov files
    if (pheno_file == None) and (cov_file == None):
        raise Exception("Sorry, you need to specify at least one of either the phenotype file or covariate file") 

    # load phenotypes
    if pheno_file == None:
        print("No separate phenotype file specified.")
        phenotypes =[]
        
    else:        
        try:
            df = pd.read_table(pheno_file, sep = " ", header = 0)
            df.columns = [col_name.lower() for col_name in df.columns]
            phenotypes = df.columns[2:]
            
        except FileNotFoundError:
            print("Specified phenotype file is not found or cannot be loaded")
            # create empty list of phenotypes 
            phenotypes = []
        
    # read in covariates if nonnull
    if cov_file == None:
        print("No covariates file specified.")
        covariates = []
        
    else: 
        try:
            cov_selected = pd.read_table(cov_file, sep = " ", header=0)
            cov_selected.columns = [col_name.lower() for col_name in cov_selected.columns]
            covariates = cov_selected.columns[2:]
            try:
                df = pd.merge(cov_selected, df, on = ["fid", "iid"])
            except :
                df = cov_selected
        except FileNotFoundError:
            print("Specified covariate file is not found or cannot be loaded")
            # create empty list of covariates to return
            covariates= []
            
    # read in pcs if nonnull
    if PC_file == None:
        print("No PC file specified.")

    else: 
        try:
            PCs = pd.read_table(PC_file, sep= " ", header=None)
            PCs.columns = ["fid", "iid"] + ["pc_" + str(s) for s in range(1, PCs.shape[1]-1)]
            df = pd.merge(df, PCs, on=["fid", "iid"])
            
        except FileNotFoundError:
            print("Specified PC file is not found or cannot be loaded")
            
    
    # return the full dataframe as well as names for covariates and phenotypes
    return df, covariates, phenotypes  


# %% Read GRM
def load_everything(prefix, pheno_file, cov_file=None, PC_file=None, k=0):
    """
    Load all covariates, phenotypes, and the GRM

    Parameters
    ----------
    prefix : string
        path to grm files.
    pheno_file : string
        path to phenotype file.
    cov_file : string, optional
        Path to covariate file. The default is None.
    PC_file : string, optional
        path to PC's file. The default is None.
    k : int, optional
        numbero f partitions for the estimate. The default is 0.

    Returns
    -------
    a tuple of the full dataframe, covariate names, phenotype names, GRM without missing vlaues, and ids.

    """
    
    print("Reading GRM: ", prefix)
    
    # Time reading the GRM and other data
    start_read = timeit.default_timer()
    
    # Read in grm
    G = ReadGRMBin(prefix)
    # Get specific detials about the GRM
    ids = G['id']
    n_phen_nona = G['n_phen_nona']
    GRM_array_nona = np.zeros((n_phen_nona, n_phen_nona))
    GRM_array_nona[np.diag_indices(n_phen_nona)] = G['diag']
    
    ###############################
    # Don't know what this is doing
    if(k == 0):
        k = n_phen_nona
    temp_i = 0
    temp = 0
    # k= args.k
    
    l = list(range(k, n_phen_nona, k))
    l.append(n_phen_nona)
    for i in l:
        cor = multirange(range(temp_i, i))
        GRM_array_nona[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
        GRM_array_nona.T[cor['b'], cor['a']] = G['off'][temp:temp+len(cor['b'])]
        temp = temp + len(cor['b'])
        del(cor)
        temp_i = i
    ################################
    
    
    df, covariates, phenotypes = load_data(pheno_file=pheno_file, cov_file=cov_file, PC_file=PC_file)
    end_read = timeit.default_timer()
    read_time = end_read - start_read
    
    print("It took " + str(read_time) + " (s) to read GRM, covariates, and phenotypes")
    print("Phenos + Covars:", df.columns)
    
    return df, covariates, phenotypes, GRM_array_nona, ids 
