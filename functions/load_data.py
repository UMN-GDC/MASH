import os
import pandas as pd
import numpy as np
from struct import unpack, calcsize
import timeit
from functions.loading_extracting_niis import load_extract_niis

#%%





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
    ids = pd.DataFrame(np.loadtxt(IDFileName, delimiter = '\t', dtype = str))
    u = np.tril_indices(ids.shape[0])
    n = len(ids.index)
    ## Read relatedness values
    grm = np.fromfile(BinFileName, dtype = dt)
    # seed empty grm
    GRM = np.zeros((n, n), dtype = dt)
    # make an upper triangle matrix
    GRM[u] = grm
    # Make the rest of the symmetric matrix
    GRM = GRM + GRM.T - np.diag(np.diag(GRM))
    return GRM


# To check if a file is missing a header
def check_header(filename):
        with open(filename) as f:
            first = f.read(1)
        return first not in '.-0123456789'


def data_loader(file) :
    # if filepath is empty it's the GRM
    if os.path.splitext(file)[-1] == "":
        df = pd.DataFrame(np.loadtxt(file+ ".grm.id", delimiter = '\t', dtype = str))
        df.columns = ["fid", "iid"]

    # check if it's the pheno or covar file
    elif check_header(file) :
        df = pd.read_table(file, sep = "\s+", header = 0)
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

def load_tables(list_of_files) :
    # load first dataframe
    df = data_loader(list_of_files[0])
    # load the rest of the data
    for file in list_of_files[1:] :
        if file != None: 
            newdf = data_loader(file)
            # merge always using the left keys such that it always aligns with the GRM
            df = pd.merge(df, newdf, on = ["fid", "iid"], how = "left")
    return df




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
    a tuple of the full dataframe, GRM without missing vlaues

    """
    
    print("Reading GRM: ", prefix)
    
    # Time reading the GRM and other data
    start_read = timeit.default_timer()
    
    # Read in grm
    GRM = ReadGRMBin(prefix)
    list_of_files = [prefix, pheno_file, PC_file, cov_file]
    df = load_tables(list_of_files)

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

