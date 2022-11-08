import os
import pandas as pd
import numpy as np
import timeit

#%%


def ReadGRMBin(prefix):
    """
    Read GCTA style binary GRM file sets into memory.

    Parameters
    ----------
    prefix : string
        filepath common to all files of the GRM that is to be read.

    Returns
    -------
    ids : pandas dataframe
        Dataframe containing the FID and IID of subjects in the same order as the GRM
    GRM : numpy array
        GRM as an (nxn) array.

    """
    print("Reading GRM: ", prefix)
    
    # Specify information about binary GRM format
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file

    # Read IDs
    ids = pd.DataFrame(np.loadtxt(prefix + ".grm.id",
                                  delimiter = '\t', dtype = str), columns = ["fid", "iid"], dtype = str)
    n = ids.shape[0]

    ## Read GRM from binary 
    grm = np.fromfile(prefix + ".grm.bin", dtype = dt)
    # seed empty grm
    GRM = np.zeros((n, n), dtype = dt)
    # make a lower triangle matrix
    l = np.tril_indices(n)
    GRM[l] = grm
    # Make the rest of the symmetric matrix
    GRM = GRM + GRM.T - np.diag(np.diag(GRM))
    return ids, GRM


# To check if a file is missing a header
def check_header(filename):
        with open(filename) as f:
            first = f.read(1)
        return first not in '.-0123456789'


def data_loader(file) :
    # check if it's the pc or covar or pheno file  (pc file won't have a header)
    if check_header(file) :
        df = pd.read_table(file, sep = "\s+", header = 0)
        # standardize column names by making them all lowercase
        df.columns = [col_name.lower() for col_name in df.columns]

    # if not it's the PC file
    else:
        df = pd.read_table(file, sep= "\s+", header=None)
        df.columns = ["fid", "iid"] + ["pc_" + str(s) for s in range(1, df.shape[1]-1)]
        df.fid = df.fid.astype("Int64")

    # make sure fid and iid are objects to join with the ids from the GRM
    df["fid"] = df.fid.astype(str)
    df["iid"] = df.iid.astype(str)
    return df

def load_tables(ids, list_of_files) :
    # load the rest of the data
    for file in list_of_files :
        if file != None: 
            newdf = data_loader(file)
            # merge always using the left keys such that it always aligns with the GRM
            ids = pd.merge(ids, newdf, on = ["fid", "iid"], how = "left")
    return ids




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

    Returns
    -------
    a tuple of the full dataframe, GRM without missing vlaues where the order of the FIDS and IIDs match between the df and GRM

    """
    
    print("Reading GRM: ", prefix)
    
    # Time reading the GRM and other data
    start_read = timeit.default_timer()
    
    # Read in grm
    ids, GRM = ReadGRMBin(prefix)
    list_of_files = [pheno_file, PC_file, cov_file]
    df = load_tables(ids, list_of_files)

    end_read = timeit.default_timer()
    read_time = end_read - start_read
    
    print("It took " + str(read_time) + " (s) to read GRM, covariates, and phenotypes")
    print("Phenos + Covars:", df.columns)
   
    # Get the phenotype names
    phenotypes = pd.read_table(pheno_file, sep = "\s+", header = 0, nrows= 0).columns.tolist()
    phenotypes = [phenotype.lower() for phenotype in phenotypes]
    phenotypes.remove("fid")
    phenotypes.remove("iid")
    return df, GRM, phenotypes 

