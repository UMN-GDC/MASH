import pandas as pd
import numpy as np
from struct import unpack, calcsize
from functions.loading_extracting_niis import load_extract_niis
#%%



def sum_n_vec(n):
    s = [int(0)] * n
    for i in range(n):
        s[i] = int(((i + 1) * (i + 2) / 2) - 1)
    return(s)


def ReadGRMBin(prefix, AllN = False):
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
    ids = ids.rename(columns={0: "FID", 1: "IID"})
    ids = ids.dropna()
    # ids["FID"] = ids.FID.astype(int)
    n_phen_nona = ids.shape[0]
    n_phen_nona = grm[i].size
    GRM_array_nona = np.zeros((n_phen_nona, n_phen_nona))
    val = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids,'N':N, "n_phen_nona" : n_phen_nona}
    return(val)


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
    return(val)




# Read data function that can load csv pheno and txt file types
# def read_datas(file_path, IDs) :
#     if(file_path.split(".")[-1] == "csv"):
#         # This is expected ot have column names
#         dat = pd.read_csv(file_path)
#         # dat.columns = ["FID", "IID"] + ["Covar_" + str(s) for s in range(1, dat.shape[1]-1)]
#     elif(file_path.split(".")[-1] == "txt"):
#         dat = pd.read_table(file_path, sep = " ")
#         # dat.columns = ["FID", "IID"] + ["Covar_" + str(s) for s in range(1, dat.shape[1] -1)]
#     elif(file_path.split(".")[-1] == "phen"):
#         dat = pd.read_table(file_path, sep = " ")
#         # dat.columns = ["FID", "IID"] + ["Pheno_" + str(s) for s in range(1, dat.shape[1]-1)]
#     elif(file_path.split(".")[-1] == "eigenvec"):
#         dat = pd.read_table(file_path, sep = " " , header=None)
#         dat.columns = ["FID", "IID"] + ["PC_" + str(s) for s in range(1, dat.shape[1] -1)]
#     elif(file_path.split(".")[-1] == "files"):
#         dat = load_extract_niis(file_path, IDs)
#         dat.columns = ["FID", "IID"] + ["Pheno_" + str(s) for s in range(1, dat.shape[1]-1 )]
#         # remove the unintentional columns that sometimes happen with phenotype and csv filetypes
#         dat = dat[dat.columns.drop(list(dat.filter(regex='Unnamed')))]
#         # dat = dat.rename(columns={0 : "FID", 1 : "IID"})
#     return(dat)

# def read_datas(file_path, IDs) :
#     pd.read_csv(file_path)


# Read covariates, PC's, and phenotype all at once
def load_data(pheno_file, cov_file=None, PC_file= None) :
    # load phenotypes
    df = pd.read_table(pheno_file, sep = " ", header = 0)
    phenotypes = df.columns
    # read in covariates if nonnull
    try:
        cov_selected = pd.read_table(cov_file, sep = " ", header=0)
        df = pd.merge(cov_selected, df, on = ["FID", "IID"])
    except:
        print("No covariates file specified or specified file is not found or cannot be loaded.")
        # onlyt load pcs if non null
    try:
        PCs = pd.read_table(PC_file, sep= " ", header=0)
        PCs.columns = ["FID", "IID"] + ["PC_" + str(s) for s in range(1, PCs.shape[1]-1)]
        df = pd.merge(PCs, df, on=["FID", "IID"])
    except:
        print("No PC file specified or specified file is not found or cannot be loaded.")
        #  if(PCs != None) :
            #    print("You  specified a PC file, without specifying how many PC's, here we assume keeping 0 PC's")
    # return the full dataframe as well as names for covariates and phenotypes
    return(df, cov_selected.columns[2:], phenotypes[2:] )


