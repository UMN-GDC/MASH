import pandas as pd
import numpy as np
from struct import unpack, calcsize



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
    val = {'diag': grm[i], 'off': np.delete(grm, i),'id': ids,'N':N}
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
def read_datas(file_path) :
 if(file_path.split(".")[-1] == "csv"):
  # read csv with write optinos
  dat = pd.read_csv(file_path, header=None)
  # Rename columns to make them neat
  n=list(range(1, dat.shape[1] -1))
  dat.columns = ["FID", "IID"] + ["Covar_" + str(s) for s in n]
 elif(file_path.split(".")[-1] == "phen"):
  dat = pd.read_table(file_path, sep = " " , header=None)
  n = list(range(1, dat.shape[1] -1))
  dat.columns = ["FID", "IID"] + ["Pheno_" + str(s) for s in n]
 elif(file_path.split(".")[-1] == "txt"):
  dat = pd.read_table(file_path, sep = " " , header=None)
  n = list(range(1, dat.shape[1] -1))
  dat.columns = ["FID", "IID"] + ["Covar_" + str(s) for s in n]
 elif(file_path.split(".")[-1] == "eigenvec"):
   dat = pd.read_table(file_path, sep = " " , header=None)
   n = list(range(1, dat.shape[1] -1))
   dat.columns = ["FID", "IID"] + ["PC_" + str(s) for s in n]

 # remove the unintentional columns that sometimes happen with phenotype and csv filetypes
 dat = dat[dat.columns.drop(list(dat.filter(regex='Unnamed')))]
 # dat = dat.rename(columns={0 : "FID", 1 : "IID"})
 return(dat)


# Read covariates, PC's, and phenotype all at once
def load_data(pheno_file, cov_file=None, PC_file=None, npc=-9) :
  # load phenotypes
  df = read_datas(pheno_file)
  # read in covariates if nonnull
  try:
    cov_selected = read_datas(cov_file)
    df = pd.merge(cov_selected, df, on = ["FID", "IID"])
  except:
    print("No covariates file specified or specified file is not found or cannot be loaded.")

  # onlyt load pcs if non null
  try:
    PCs = read_datas(PC_file)
    PCs= PCs.iloc[:, list(range(npc + 2))]
    df = pd.merge(PCs, df, on=["FID", "IID"])
    
  except:
    print("No PC file specified or specified file is not found or cannot be loaded.")
  if(npc == -9 and PC != None) :
    print("You  specified a PC file, without specifying how many PC's, here we assume keeping 0 PC's")
  return(df)


