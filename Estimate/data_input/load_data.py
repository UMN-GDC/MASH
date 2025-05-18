import os
import pandas as pd
import numpy as np
import timeit

#%%


def ReadGRMBin(prefix, sub_ids = None):
    """
    Read GCTA style binary GRM file sets into memory.

    Parameters
    ----------
    prefix : string
        filepath common to all files of the GRM that is to be read.
    sub_ids : str
        OPTIONAL: filepath to FIDs/IIDS of subset of the sample
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
    ids = pd.read_table(prefix + ".grm.id", names = ["FID", "IID"], dtype = str)
    ids["missing"] = ids["FID"].isna() | ids["IID"].isna()
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

    # drop missing
    GRM = GRM[np.invert(ids["missing"]),:][:,np.invert(ids["missing"])]
    ids = ids.dropna()[["FID", "IID"]]
    ids.FID = ids.FID.astype(int)
    
    if sub_ids != None :
        ids2 = pd.read_table(sub_ids, names= ["FID", "IID"], dtype = str)
        # keep the ids that overlap with the additionally specified ids
        ids = ids.reset_index().merge(ids2, on = ["FID", "IID"]).set_index("index")
        # Subset the GRM too
        GRM = GRM[ids.index, :][:, ids.index]        

    return ids, GRM

def insert_underscore(s):
    return s[:4] + '_' + s[4:]


def load_tables(ids= None, args = None) :
    # load the rest of the data
    if args["covar"] != None :
        covDF = pd.read_table(args["covar"], sep = "\s+")
    if "ids" in locals() :
        ids = covDF[["FID", "IID"]] 
    if args["PC"] != None :
        pcDF = pd.read_table(args["PC"], sep = "\s+", header= None)
        pcDF.columns = ["FID", "IID"] + ["pc_" + str(s) for s in range(1, pcDF.shape[1]-1)]
        pcDF = pcDF.dropna()
        pcDF.FID = pcDF.FID.astype(int)

    
    if os.path.splitext(args["pheno"])[1] == ".parquet":

        pheno=pd.read_parquet(args["pheno"]).reset_index(names="IID")
        pheno.IID = pheno.IID.astype(str)
        pheno['IID'] = pheno['IID'].apply(insert_underscore)               
    else :
        pheno=pd.read_table(args["pheno"], sep = "\s+")

    for dframe in ["covDF", "pcDF", "pheno"] :
        if dframe in locals() :
            ids = pd.merge(ids, eval(dframe), on = ["FID", "IID"], how = "left")

    return ids




# %% Read GRM
def load_everything(args, k=0):
    """
    Load all covariates, phenotypes, and the GRM

    Parameters
    ----------
    args : dict
        path to grm files.

    Returns
    -------
    a tuple of the full dataframe, GRM without missing vlaues where the order of the FIDS and IIDs match between the df and GRM

    """
    
    print("Reading GRM: ", args["prefix"])
    
    # Time reading the GRM and other data
    start_read = timeit.default_timer()
    
    # Read in grm
    ids, GRM = ReadGRMBin(args["prefix"], args["ids"])
    df = load_tables(ids, args)

    end_read = timeit.default_timer()
    read_time = end_read - start_read
    
    print("It took " + str(read_time) + " (s) to read GRM, covariates")
    print("Phenos + Covars:", df.columns.tolist())
   
    # Get the phenotype names
    if os.path.splitext(args["pheno"])[1] == ".parquet":

        phenotypes=pd.read_parquet(args["pheno"]).T.rename(columns=lambda x: 'o' + str(x)).reset_index(names="IID").columns.tolist()
        phenotypes.remove("IID")
    else: 
        phenotypes = pd.read_table(args["pheno"], sep = "\s+", header = 0, nrows= 0).columns.tolist()
        phenotypes.remove("FID")
        phenotypes.remove("IID")
    print(df.shape)
    ids = ids.dropna()
    ids["FID"] = ids.FID.astype(int)
    return df, GRM, phenotypes, ids

