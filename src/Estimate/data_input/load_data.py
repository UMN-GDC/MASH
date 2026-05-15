import os
import logging
import pandas as pd
import numpy as np
import timeit

logger = logging.getLogger(__name__)

#%%


def ReadGRMBin(prefix, sub_ids = None, args = None):
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
    iid_col = args.get("iid_col", "IID") if args else "IID"
    fid_col = args.get("fid_col", "FID") if args else "FID"

    print("Reading GRM: ", prefix)

    # Specify information about binary GRM format
    dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file

    # Read IDs - handle custom column names
    ids = pd.read_table(prefix + ".grm.id", header=None, dtype = str)
    ids.columns = ids.columns.map({0: fid_col, 1: iid_col})
    ids = rename_id_cols(ids, args) if args else ids
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

    if sub_ids != None :
        sub_fid_col = args.get("fid_col", "FID") if args else "FID"
        sub_iid_col = args.get("iid_col", "IID") if args else "IID"
        ids2 = pd.read_table(sub_ids, header=None, dtype = str)
        ids2.columns = ids2.columns.map({0: sub_fid_col, 1: sub_iid_col})
        ids2 = rename_id_cols(ids2, args) if args else ids2
        # keep the ids that overlap with the additionally specified ids
        ids = ids.reset_index().merge(ids2, on = ["FID", "IID"]).set_index("index")
        # Subset the GRM too
        GRM = GRM[ids.index, :][:, ids.index]

    return ids, GRM

def insert_underscore(s):
    return s[:4] + '_' + s[4:]


def find_col(header, col_name, desc=""):
    idx = [i for i, c in enumerate(header) if str(c).upper() == str(col_name).upper()]
    if len(idx) == 0:
        raise ValueError(f"Column '{col_name}' not found in {desc}. Available: {header}")
    return idx[0]


def rename_id_cols(df, args):
    """Rename FID/IID columns to canonical names based on args."""
    if df is None:
        return df
    fid_col = args.get("fid_col", "FID")
    iid_col = args.get("iid_col", "IID")
    renames = {}
    for old, new in [(fid_col, "FID"), (iid_col, "IID")]:
        if old != new and old in df.columns:
            renames[old] = new
    if renames:
        df = df.rename(columns=renames)
    return df


def validate_input_files(args):
    """
    Validate that all required input files exist and are readable.
    
    Parameters
    ----------
    args : dict
        Dictionary of arguments containing file paths
        
    Raises
    ------
    FileNotFoundError
        If any required file doesn't exist
    """
    required_files = []
    optional_files = []
    
    # Check prefix (required for most methods)
    if args.get("prefix") and args["prefix"] != "None":
        prefix = args["prefix"]
        required_files.append((f"{prefix}.grm.bin", "GRM binary"))
        required_files.append((f"{prefix}.grm.id", "GRM IDs"))
        # .grm.N.bin is optional
    
    # Check phenotype file (required)
    pheno_input = args.get("pheno")
    if pheno_input is not None:
        if isinstance(pheno_input, list):
            for pf in pheno_input:
                required_files.append((pf, "Phenotype"))
        elif pheno_input != "None":
            required_files.append((pheno_input, "Phenotype"))
    
    # Check optional files
    if args.get("PC") and args["PC"] != "None":
        pc_input = args["PC"]
        if isinstance(pc_input, list):
            for pc_file in pc_input:
                optional_files.append((pc_file, "PC"))
        else:
            optional_files.append((pc_input, "PC"))
    
    if args.get("covar") and args["covar"] != "None":
        covar_input = args["covar"]
        if isinstance(covar_input, str):
            covar_files = [f.strip() for f in covar_input.split(',')]
        else:
            covar_files = list(covar_input)
        for cf in covar_files:
            optional_files.append((cf, "Covariate"))
    
    # Validate required files
    missing = []
    for filepath, filetype in required_files:
        if not os.path.isfile(filepath):
            missing.append((filepath, filetype))
    
    if missing:
        for filepath, filetype in missing:
            logger.error(f"Required {filetype} file not found: {filepath}")
        raise FileNotFoundError(f"Missing required files: {[f[0] for f in missing]}")
    
    # Log optional files
    for filepath, filetype in optional_files:
        if os.path.isfile(filepath):
            logger.info(f"Using {filetype} file: {filepath}")
        else:
            logger.warning(f"Optional {filetype} file not found: {filepath}")
    
    # Validate GRM files
    if args.get("prefix") and args["prefix"] != "None":
        prefix = args["prefix"]
        if not os.path.isfile(f"{prefix}.grm.bin"):
            raise FileNotFoundError(f"GRM binary file not found: {prefix}.grm.bin")
        if not os.path.isfile(f"{prefix}.grm.id"):
            raise FileNotFoundError(f"GRM ID file not found: {prefix}.grm.id")
    
    logger.info("All input files validated successfully")


def validate_column_types(df, filepath, file_type):
    """
    Validate and convert column types for key data columns.
    
    Parameters
    ----------
    df : pandas DataFrame
        DataFrame to validate
    filepath : str
        Path to the file (for logging)
    file_type : str
        Type of file (PC, phenotype, covariate)
        
    Returns
    -------
    pandas DataFrame with validated types
    """
    if file_type == "PC":
        # PCs should be numeric
        pc_cols = [c for c in df.columns if c.lower().startswith("pc")]
        logger.debug(f"PC file columns: {df.columns.tolist()}")
        logger.debug(f"PC cols found: {pc_cols}")
        for col in pc_cols:
            if col not in ["FID", "IID"]:
                try:
                    col_data = df[col]
                    if not isinstance(col_data, pd.Series):
                        logger.warning(f"PC column '{col}' in {filepath} is not a Series (type: {type(col_data)}), skipping numeric conversion")
                        continue
                    df[col] = pd.to_numeric(col_data, errors='raise')
                except Exception as e:
                    logger.warning(f"PC column '{col}' in {filepath} contains non-numeric values, converting to NaN: {e}")
                    df[col] = pd.to_numeric(col_data, errors='coerce')
        if pc_cols:
            logger.debug(f"Validated {len(pc_cols)} PC columns as numeric")
    
    elif file_type == "phenotype":
        # Phenotype columns should be numeric
        pheno_cols = [c for c in df.columns if c not in ["FID", "IID"]]
        for col in pheno_cols:
            if col not in ["FID", "IID"]:
                try:
                    df[col] = pd.to_numeric(df[col], errors='raise')
                except Exception as e:
                    logger.warning(f"Phenotype column '{col}' in {filepath} contains non-numeric values, converting to NaN: {e}")
                    df[col] = pd.to_numeric(df[col], errors='coerce')
        if pheno_cols:
            logger.debug(f"Validated {len(pheno_cols)} phenotype columns as numeric")
    
    elif file_type == "covariate":
        # For covariates, try to convert to numeric but don't warn for string columns (they may be categorical)
        for col in df.columns:
            if col not in ["FID", "IID"]:
                # Try to convert, keep as-is if it fails (likely categorical)
                try:
                    df[col] = pd.to_numeric(df[col], errors='raise')
                except (ValueError, TypeError):
                    pass  # Keep as string - likely categorical
    
    return df


def _read_delimited_file(filepath, has_header=True, default_cols=None, args=None):
    """
    Read a file with unambiguous delimiter based on extension.

    Returns DataFrame with FID, IID columns.
    """
    import os
    fid_col = args.get("fid_col", "FID") if args else "FID"
    iid_col = args.get("iid_col", "IID") if args else "IID"
    ext = os.path.splitext(filepath)[1].lower()

    # Determine separator based on extension
    if ext in ['.tsv', '.tab']:
        sep = '\t'
    elif ext == '.csv':
        sep = ','
    elif ext in ['.phe', '.pheno', '.phen', '.txt', '.covar']:
        sep = None  # whitespace
    else:
        # Default: try whitespace (plink format)
        sep = None
    
    # Handle case where sep is None but we need to pass a valid separator to pandas
    if sep is None:
        sep = r'\s+'  # whitespace regex for pandas
        
    # Log separator for debugging
    logger.debug(f"Reading file {filepath} with sep={repr(sep)} (type: {type(sep)})")
    
    # Validate sep is not None before passing to pandas
    if sep is None:
        logger.error("sep is None after processing! Setting to default whitespace.")
        sep = r'\s+'
        
    # Additional validation - make sure sep is actually a string
    if not isinstance(sep, str):
        logger.error(f"sep is not a string! Got {type(sep)}: {sep}. Setting to default whitespace.")
        sep = r'\s+'

    logger.debug(f"Final sep value: {repr(sep)} (type: {type(sep)})")
    
    # Try to read with header
    try:
        logger.debug(f"About to call pd.read_table with filepath={repr(filepath)}, sep={repr(sep)}, header={None if not has_header else 0}")
        # Ensure sep is a valid string for pandas
        if sep is None:
            logger.error("sep is None! Forcing to whitespace separator.")
            sep = r'\s+'
        elif not isinstance(sep, str):
            logger.error(f"sep is not a string! Got {type(sep)}: {sep}. Forcing to whitespace separator.")
            sep = r'\s+'
        logger.debug(f"Final sep value for pd.read_table: {repr(sep)}")
        df = pd.read_table(filepath, sep=sep, header=None if not has_header else 0)
    except Exception as e:
        logger.error(f"Failed to read {filepath} with sep={repr(sep)}: {e}")
        logger.error(f"sep type: {type(sep)}, sep value: {sep}")
        # Log more details about the exception
        import traceback
        logger.error(f"Traceback: {traceback.format_exc()}")
        # Try with default sep=None as fallback
        logger.info(f"Trying fallback with sep=None")
        df = pd.read_table(filepath, sep=None, header=None if not has_header else 0)

    # Strip whitespace from column names (common issue with some files)
    if has_header:
        df.columns = df.columns.str.strip()
    
    # Normalize #FID/#IID hashed column names in the first two columns
    if has_header:
        renames = {}
        for col in list(df.columns):
            stripped = col.lstrip('#')
            if stripped != col and stripped.upper() in ['FID', 'IID']:
                renames[col] = stripped
        if renames:
            df = df.rename(columns=renames)

    # Check if FID/IID columns exist with custom names (case-insensitive match)
    def has_col(name):
        return any(str(c).upper() == str(name).upper() for c in df.columns)

    if has_col(fid_col) and has_col(iid_col):
        # Ensure FID and IID are string type for consistent merging
        df[fid_col] = df[fid_col].astype(str)
        df[iid_col] = df[iid_col].astype(str)
        df = rename_id_cols(df, args)
        return df

    # Fallback to canonical FID/IID columns (in case config uses custom names but file has standard names)
    if has_col('FID') and has_col('IID'):
        df['FID'] = df['FID'].astype(str)
        df['IID'] = df['IID'].astype(str)
        df = rename_id_cols(df, args)
        return df

    # Handle alternative column names: participant_id -> IID (for files without FID)
    # Keep session_id as a separate column if present
    if has_col('participant_id') and not has_col(fid_col):
        col_rename = {'participant_id': iid_col}
        df = df.rename(columns=col_rename)
        df[fid_col] = df[iid_col].astype(str)
        df[iid_col] = df[iid_col].astype(str)
        df = rename_id_cols(df, args)
        return df

    # Handle subjectkey -> IID (common in some phenotype files)
    if has_col('subjectkey') and not has_col(fid_col):
        col_rename = {'subjectkey': iid_col}
        df = df.rename(columns=col_rename)
        df[fid_col] = df[iid_col].astype(str)
        df[iid_col] = df[iid_col].astype(str)
        df = rename_id_cols(df, args)
        return df

    # FID/IID were not found or no header provided - column names were inferred
    logger.info(f"File {filepath}: FID/IID columns not found, inferring from first two columns")

    # No header or doesn't have FID/IID - assign columns
    n_cols = df.shape[1]
    if default_cols is not None and len(default_cols) > 0:
        # Use default_cols, fill in with generated names if needed
        col_names = list(default_cols)
        extra_cols_needed = n_cols - len(default_cols)
        if extra_cols_needed > 0:
            # Generate additional column names
            base_name = "col_"
            if default_cols:
                last_col = default_cols[-1]
                if last_col.startswith("pheno_"):
                    base_name = "pheno_"
                    try:
                        last_num = int(last_col.split("_")[-1])
                        col_names = default_cols + [f"pheno_{i}" for i in range(last_num + 1, last_num + extra_cols_needed + 1)]
                    except:
                        pass
                elif last_col.startswith("pc_"):
                    base_name = "pc_"
                elif last_col.startswith("cov_"):
                    base_name = "cov_"
            if len(col_names) < n_cols:
                col_names = default_cols + [f"{base_name}{i}" for i in range(1, extra_cols_needed + 1)]
        col_names = col_names[:n_cols]
    else:
        col_names = ["FID", "IID"] + [f"col_{i}" for i in range(n_cols - 2)]
    df.columns = col_names

    # Ensure FID and IID are string type for consistent merging
    if 'FID' in df.columns:
        df['FID'] = df['FID'].astype(str)
    if 'IID' in df.columns:
        df['IID'] = df['IID'].astype(str)

    return df


def _apply_filter(df, filter_expr):
    """
    Apply a filter expression to a dataframe.
    
    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to filter
    filter_expr : str
        Filter expression like "age>30", "sex==1", "site!=2"
    
    Returns
    -------
    pandas.DataFrame
        Filtered dataframe
    """
    import re
    
    if df is None or filter_expr is None:
        return df
    
    # Supported operators
    operators = ['==', '!=', '>=', '<=', '>', '<']
    op = None
    for o in operators:
        if o in filter_expr:
            op = o
            break
    
    if op is None:
        logger.warning(f"Invalid filter expression: {filter_expr}. No operator found.")
        return df
    
    parts = filter_expr.split(op)
    if len(parts) != 2:
        logger.warning(f"Invalid filter expression: {filter_expr}")
        return df
    
    col_name = parts[0].strip()
    value_str = parts[1].strip()
    
    if col_name not in df.columns:
        logger.warning(f"Column '{col_name}' not found in data. Available: {df.columns.tolist()}")
        return df
    
    # Try to parse the value (numeric or string)
    try:
        # Try numeric first
        value = float(value_str)
        # If it's an integer, keep it as is for comparison
        if value == int(value):
            value = int(value)
    except ValueError:
        # Keep as string (strip quotes if present)
        value = value_str.strip("'\"")
    
    # Apply filter
    original_count = len(df)
    if op == '==':
        df = df[df[col_name] == value]
    elif op == '!=':
        df = df[df[col_name] != value]
    elif op == '>':
        df = df[df[col_name] > value]
    elif op == '<':
        df = df[df[col_name] < value]
    elif op == '>=':
        df = df[df[col_name] >= value]
    elif op == '<=':
        df = df[df[col_name] <= value]
    
    filtered_count = original_count - len(df)
    if filtered_count > 0:
        logger.info(f"Filtered {filtered_count} rows from {original_count} ({filter_expr})")
    else:
        logger.info(f"Filter applied: {filter_expr}")
    
    return df


def load_tables(ids= None, args = None) :
    # load the rest of the data
    fid_col = args.get("fid_col", "FID") if args else "FID"
    iid_col = args.get("iid_col", "IID") if args else "IID"
    
    # Handle multiple covariate files (list or comma-separated string)
    if args.get("covar") is not None:
        covar_input = args["covar"]
        logger.debug(f"Processing covar input: {covar_input} (type: {type(covar_input)})")
        if isinstance(covar_input, str):
            covar_files = [f.strip() for f in covar_input.split(',')]
            logger.debug(f"Split covar string into files: {covar_files}")
        elif isinstance(covar_input, list):
            covar_files = covar_input
            logger.debug(f"Using covar list as files: {covar_files}")
        else:
            covar_files = [covar_input]
            logger.debug(f"Wrapped single covar input in list: {covar_files}")
        
        covDFs = []
        for cov_file in covar_files:
            logger.debug(f"Processing covariate file: {cov_file}")
            try:
                cov_df = _read_delimited_file(cov_file, default_cols=["FID", "IID", "cov_1"], args=args)
                logger.debug(f"Successfully read {cov_file}, shape: {cov_df.shape}")
                # Drop session_id if present (used for filtering, not needed as covariate)
                if 'session_id' in cov_df.columns:
                    cov_df = cov_df.drop(columns=['session_id'])
                # Apply covariate filter to each file BEFORE merging to avoid duplicates
                # Only apply if the filter column exists in this file
                if args.get("covar_filter"):
                    filter_expr = args["covar_filter"]
                    logger.debug(f"Applying filter: {filter_expr}")
                    # Extract column name from filter expression
                    for op in ['==', '!=', '>=', '<=', '>', '<']:
                        if op in filter_expr:
                            filter_col = filter_expr.split(op)[0].strip()
                            break
                    logger.debug(f"Filter column: {filter_col}")
                    if filter_col in cov_df.columns:
                        logger.debug(f"Filter column {filter_col} found in dataframe")
                        cov_df = _apply_filter(cov_df, filter_expr)
                        logger.debug(f"After filtering, shape: {cov_df.shape}")
                    else:
                        logger.warning(f"Filter column '{filter_col}' not found in {cov_file}, skipping filter for this file")
                covDFs.append(cov_df)
                logger.debug(f"Added {cov_file} to covDFs list, now {len(covDFs)} files processed")
            except Exception as e:
                logger.warning(f"Could not read covariate file {cov_file}: {e}")
                import traceback
                logger.warning(f"Traceback: {traceback.format_exc()}")
        
        if covDFs:
            # Merge all covariate dataframes, keeping unique columns
            covDF = covDFs[0]
            for df in covDFs[1:]:
                # Find overlapping columns (excluding FID, IID)
                overlap = set(covDF.columns) & set(df.columns) - {"FID", "IID"}
                if overlap:
                    # Rename overlapping columns in the second df before merging
                    df = df.rename(columns={c: f"{c}_2" for c in overlap})
                covDF = pd.merge(covDF, df, on=["FID", "IID"], how="inner")
    else:
        covDF = None
    
    # Note: covar_filter already applied above per-file to avoid duplicates
    
    if args.get("PC") != None :
        # Auto-detect if the PC file has a header by checking the first line
        try:
            with open(args["PC"], 'r') as f:
                first_line = f.readline().strip()
            first_tokens = first_line.split()
            # Heuristic: has a header if first two tokens are FID/IID-like and
            # at least one later token starts with PC
            has_pc_header = (
                len(first_tokens) >= 3 and
                first_tokens[0].upper() in ['FID', '#FID', 'FAM', '#FAM'] and
                first_tokens[1].upper() in ['IID', '#IID', 'ID', '#ID'] and
                any(t.upper().startswith('PC') for t in first_tokens[2:])
            )
        except Exception:
            has_pc_header = False

        if has_pc_header:
            # File has a proper header — read normally and lowercase PC columns
            pcDF = _read_delimited_file(args["PC"], args=args)
            # Handle #FID prefix if present
            if pcDF.columns[0].startswith('#'):
                pcDF = pcDF.rename(columns={'#FID': 'FID'})
            # Lowercase PC column names
            col_mapping = {c: c.lower() for c in pcDF.columns if c.upper().startswith('PC') and c not in ['FID', 'IID']}
            if col_mapping:
                pcDF = pcDF.rename(columns=col_mapping)
                logger.info(f"Normalized PC column names: {col_mapping}")
            pcDF = pcDF.dropna()
        else:
            # No header — use the existing fallback logic
            pcDF = _read_delimited_file(args["PC"], default_cols=["FID", "IID", "PC1", "PC2"], has_header=False, args=args)
            if pcDF.columns[0].startswith('#'):
                pcDF = pcDF.rename(columns={'#FID': 'FID'})
            # Check if first column contains space-separated values (no header case)
            # If so, split them properly
            if 'FID' not in pcDF.columns or isinstance(pcDF.columns[0], int):
                # Need to split the first column into multiple columns
                # Format: FID IID PC1 PC2 ... PCn
                first_col = pcDF.iloc[:, 0].astype(str)
                # Split by whitespace
                split_data = first_col.str.split(expand=True)
                pcDF = pd.concat([split_data.iloc[:, 0:2], split_data.iloc[:, 2:]], axis=1)
                pcDF.columns = ['FID', 'IID'] + [f'pc{i}' for i in range(1, pcDF.shape[1]-1)]
            pcDF = pcDF.dropna()
            # Normalize PC column names - handle both "PC1" style and numeric column names
            # Rename numeric columns (3, 4, 5, ...) to pc1, pc2, pc3, ...
            new_cols = {}
            pc_count = 0
            for col in pcDF.columns:
                if col not in ['FID', 'IID'] and not col.lower().startswith('pc'):
                    pc_count += 1
                    new_cols[col] = f'pc{pc_count}'
            if new_cols:
                pcDF = pcDF.rename(columns=new_cols)
                logger.info(f"Renamed PC columns: {new_cols}")
            # Also lowercase any existing PC columns
            col_mapping = {c: c.lower() for c in pcDF.columns if c.lower().startswith("pc") and c != c.lower()}
            if col_mapping:
                pcDF = pcDF.rename(columns=col_mapping)
                logger.info(f"Normalized PC column names: {col_mapping}")
            logger.debug(f"PC dataframe columns after processing: {pcDF.columns.tolist()}")
        # Validate PC columns are numeric
        pcDF = validate_column_types(pcDF, args["PC"], "PC")
    
    # Handle phenotype - always treat as list for consistency
    pheno_input = args.get("pheno")
    if pheno_input is not None:
        # Normalize to list of pheno files
        if isinstance(pheno_input, str):
            pheno_files = [f.strip() for f in pheno_input.split(',')] if ',' in pheno_input else [pheno_input]
        elif isinstance(pheno_input, list):
            pheno_files = pheno_input
        else:
            pheno_files = [pheno_input]
        
        pheno_dfs = []
        for pf in pheno_files:
            if pf and pf != "None":
                print(f"Reading phenotype file: {pf}")
                ext = os.path.splitext(pf)[1].lower()
                
                if ext == ".parquet":
                    pdf = pd.read_parquet(pf).reset_index(names="IID")
                    pdf.IID = pdf.IID.astype(str)
                    pdf['IID'] = pdf['IID'].apply(insert_underscore)
                else:
                    pdf = _read_delimited_file(pf, default_cols=["FID", "IID"], args=args)
                
                # Drop session_id if present (used for filtering, not needed in merge)
                if 'session_id' in pdf.columns:
                    pdf = pdf.drop(columns=['session_id'])
                
                # Apply phenotype filter BEFORE validation
                if args.get("pheno_filter"):
                    pdf = _apply_filter(pdf, args.get("pheno_filter"))
                
                pdf = validate_column_types(pdf, pf, "phenotype")
                pheno_dfs.append(pdf)
        
        if pheno_dfs:
            # Merge all phenotype dataframes, keeping unique columns
            pheno = pheno_dfs[0]
            for pdf in pheno_dfs[1:]:
                overlap = set(pheno.columns) & set(pdf.columns) - {"FID", "IID"}
                if overlap:
                    pdf = pdf.rename(columns={c: f"{c}_2" for c in overlap})
                pheno = pd.merge(pheno, pdf, on=["FID", "IID"], how="inner")
            
            # After validation, drop rows where key phenotype columns are all NaN
            pheno_cols = [c for c in pheno.columns if c not in ["FID", "IID"]]
            pheno = pheno.dropna(subset=pheno_cols, how='all')
    
    # Validate covariate columns
    if covDF is not None:
        covDF = validate_column_types(covDF, args.get("covar", "unknown"), "covariate")

    for dframe in ["covDF", "pcDF", "pheno"] :
        if dframe in locals() and eval(dframe) is not None :
            df = eval(dframe).copy()
            # Deduplicate by IID to prevent row multiplication during merge
            # Keep the first occurrence (or drop duplicates)
            df = df.drop_duplicates(subset=['IID'], keep='first')
            if 'FID' not in df.columns:
                df['FID'] = df['IID']
            # Ensure consistent types
            df['FID'] = df['FID'].astype(str)
            df['IID'] = df['IID'].astype(str)
            ids['FID'] = ids['FID'].astype(str)
            # Drop FID from df to avoid duplicate column issues on merge
            if 'FID' in df.columns and 'FID' in ids.columns:
                df = df.drop(columns=['FID'])
            # Merge on IID only since FID may differ between files (e.g., GRM has 
            # numeric FID but phenotype has participant_id style FID)
            ids = pd.merge(ids, df, on = "IID", how = "left")
            # Ensure FID is preserved from GRM after merge
            if 'FID' not in ids.columns:
                ids['FID'] = ids['IID']
    
    # Drop any rows where FID or IID are NaN
    original_count = len(ids)
    drop_cols = [c for c in ["FID", "IID"] if c in ids.columns]
    ids = ids.dropna(subset=drop_cols)
    dropped_count = original_count - len(ids)
    if dropped_count > 0:
        logger.info(f"Dropped {dropped_count} rows with missing FID/IID")
    
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
    
    validate_input_files(args)
    
    print("Reading GRM: ", args["prefix"])
    
    # Time reading the GRM and other data
    start_read = timeit.default_timer()
    
    # Read in grm
    ids, GRM = ReadGRMBin(args["prefix"], args.get("ids"), args=args)
    df = load_tables(ids, args)

    end_read = timeit.default_timer()
    read_time = end_read - start_read
    
    print("It took " + str(read_time) + " (s) to read GRM, covariates")
    print("Phenos + Covars:", df.columns.tolist())
   
    # Get the phenotype names - handle list of files
    pheno_input = args.get("pheno")
    phenotypes = []
    if pheno_input is not None:
        if isinstance(pheno_input, list):
            for pf in pheno_input:
                if pf and pf != "None":
                    ext = os.path.splitext(pf)[1].lower()
                    if ext == ".parquet":
                        pdf_phenos = pd.read_parquet(pf).reset_index(names="IID").columns.tolist()
                        pdf_phenos.remove("IID")
                    else:
                        pdf = _read_delimited_file(pf, default_cols=["FID", "IID"], args=args)
                        # Drop session_id if present (same as load_tables)
                        if 'session_id' in pdf.columns:
                            pdf = pdf.drop(columns=['session_id'])
                        pdf_phenos = [c for c in pdf.columns if c not in ["FID", "IID"]]
                    phenotypes.extend(pdf_phenos)
        else:
            ext = os.path.splitext(pheno_input)[1].lower()
            if ext == ".parquet":
                phenotypes = pd.read_parquet(pheno_input).reset_index(names="IID").columns.tolist()
                phenotypes.remove("IID")
            else:
                pheno_df = _read_delimited_file(pheno_input, default_cols=["FID", "IID"], args=args)
                if 'session_id' in pheno_df.columns:
                    pheno_df = pheno_df.drop(columns=['session_id'])
                phenotypes = [c for c in pheno_df.columns if c not in ["FID", "IID"]]
    
    print(df.shape)
    ids = ids.dropna()
    return df, GRM, phenotypes, ids

